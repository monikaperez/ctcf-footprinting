'''Run mokapot.'''

import os
import argparse
import pandas as pd

import mokapot
from xgboost import XGBClassifier
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt

def main(args):

    # ensure input file exists
    pin_file = os.path.abspath(args.input_file)
    output_dir = os.path.abspath(args.output_dir)
    assert os.path.exists(pin_file), f"ERROR: input file not found: {pin_file}"
    assert os.path.exists(output_dir), f"ERROR: output directory not found: {output_dir}"

    # read in args
    save_fig = args.plot
    file_root = args.file_root
    file_root = '_'.join(filter(None, ("CTCF", file_root)))

    print("Pin feature data file: {}".format(os.path.basename(pin_file)))
    print("Saving output to: {}".format(output_dir))
    print("Using file_root: {}".format(file_root))
    print("\n")

    # load pin file
    pin = pd.read_csv(pin_file, sep="\t")

    # run mokapot
    train_fdr = 0.10
    test_fdr = 0.05
    print("Using train FDR = {}".format(train_fdr))
    print("Using test FDR = {}".format(test_fdr))
    moka_conf, models, train_psms = make_moka_model(pin, file_root, train_fdr=train_fdr, test_fdr=test_fdr, save_as=output_dir)

    if save_fig:
        plot_fdr = 0.05
        print("Using plot eval_fdr = {}".format(plot_fdr))
        # plot output
        plot_moka_res(moka_conf=moka_conf, train_psms=train_psms, save_as=output_dir, file_root=file_root, eval_fdr = plot_fdr)
    else:
        print("Skipping plots.")

    print("Woohoo! Done!")


def make_moka_model(pin, file_root, train_fdr=0.10, test_fdr=0.05, subset_max_train=2_000_000, save_as=True):
    '''Make accessibility model with mokapot.'''

    print("Making accessibility model.")
    print("dataset size: {:,}".format(pin.shape[0]))

    train = pin
    # PSM (peptide-spectrum matches, proteomics thing)
    train_psms = mokapot.read_pin(train)
    print("Train - rows: {:,} | columns: {:,}".format(train.shape[0], train.shape[1]))
    # weighing for imbalanced data
    scale_pos_weight = sum(train.Label == -1) / sum(train.Label == 1)

    # hyper paramaters for XGBoost
    grid = {
        # decision trees in model
        "n_estimators": [25, 50, 100],
        "scale_pos_weight": [scale_pos_weight],
        # depth of tree
        "max_depth": [3, 6, 9],
        # rows a terminal leaf need to touch to be kept
        "min_child_weight": [3, 6, 9],
        # learning rate
        "gamma": [0.1, 1, 10],
    }

    # set up model (ensemble classifier)
    xgb_mod = GridSearchCV(
        XGBClassifier(eval_metric = "auc"),
        param_grid = grid,
        cv = 3,
        # receiver operator curve
        scoring = "roc_auc",
        verbose = 2,
    )

    # machine learning model to re-score PSMs
    print("Making model")
    mod = mokapot.Model(xgb_mod, train_fdr = train_fdr, subset_max_train = subset_max_train)

    # run mokapot
    print("Brewing ...")
    moka_conf, models = mokapot.brew(train_psms, mod, test_fdr = test_fdr)

    # save model to file
    #model_file = "{}/CTCF_{}.mokapot.model".format(save_as, file_root)
    #models.save_model(out_file=model_file)

    if save_as is not False:
        # save results to txt file
        result_files = moka_conf.to_txt(dest_dir=save_as, file_root=file_root)
        print("Saved formatted features to: {}\n{}".format(result_files[0], result_files[1]))
    
    return (moka_conf, models, train_psms)



def plot_moka_res(moka_conf, train_psms, save_as, file_root, eval_fdr=0.05):
    '''Plot mokaplot output.'''
    
    print("Plotting mokapot results.")
    # PSM (peptide-spectrum matches, proteomics thing)
    # confidence estimates based on the best original feature
    tide_conf = train_psms.assign_confidence()

    # PSMs
    moka_psms = (moka_conf.psms["mokapot q-value"] <= eval_fdr).sum()
    tide_psms = (tide_conf.psms["mokapot q-value"] <= eval_fdr).sum()
    print(f"PSMs gained by mokapot: {moka_psms - tide_psms:,d}")

    # Peptides
    moka_peps = (moka_conf.peptides["mokapot q-value"] <= eval_fdr).sum()
    tide_peps = (tide_conf.peptides["mokapot q-value"] <= eval_fdr).sum()
    print(f"Peptides gained by mokapot: {moka_peps - tide_peps:,d}")

    moka_features_gained = [moka_psms, moka_peps]

    title = file_root

    # plot data
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    colors = ("#343131", "#24B8A0")

    # plot performance
    for ax, level in zip(axs, tide_conf.levels):
        tide_conf.plot_qvalues(level=level, c=colors[0], ax=ax,
                               label="Tide combined p-value")
        moka_conf.plot_qvalues(level=level, c=colors[1], ax=ax,
                               label="mokapot")
        ax.legend(frameon=False)

    fig.suptitle("{} - cumulative number of discoveries over q-values\n({:,} peptides gained)".format(title, moka_peps))

    plt.tight_layout()

    # save fig
    fig_file = "{}/{}.mokapot_performance.pdf".format(save_as, title,)
    print("Saving fig to: {}".format(fig_file))
    plt.savefig(fname=fig_file, bbox_inches="tight", transparent=True)


if __name__ == "__main__":
    
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True, help="Pin file with features.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output file location.")
    parser.add_argument("-r", "--file_root", required=False, default=None, 
                        help="HIGHLY recommended optional prefix for mokapot output files. This will always be CTCF_(file_root).mokapot.psms.txt.\
                              This prefix will also be used in figures.")
    parser.add_argument("-plt", "--plot", required=False, default=True, help="Save plots to output file.")
    args = parser.parse_args()
    
    # run script
    main(args)