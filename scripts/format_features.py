'''Format features to mokapot pin.'''

import os
import argparse
import pandas as pd

def main(args):

    # ensure input files exists
    pos_file = os.path.abspath(args.pos_features)
    neg_file = os.path.abspath(args.neg_features)
    output_dir = os.path.abspath(args.output_dir)
    assert os.path.exists(pos_file), f"ERROR: input file not found: {pos_file}"
    assert os.path.exists(neg_file), f"ERROR: negative file not found: {neg_file}"
    assert os.path.exists(output_dir), f"ERROR: output directory not found: {output_dir}"

    print("Positive data file: {}".format(os.path.basename(pos_file)))
    print("Negative data file: {}".format(os.path.basename(neg_file)))
    print("Saving output to: {}".format(output_dir))
    print("\n")

    # read in args
    file_root = args.file_root
    save_features = args.save_features

    # format positive features
    pin_file = pos_file.replace(".txt", ".pin") if save_features == True else None
    df_pos = format_to_pin(pos_file, "positive", save_as=pin_file)
    print("\n")
    
    # format negative features
    pin_file = neg_file.split(".")[0] + ".pin"
    df_neg = format_to_pin(neg_file, "negative", save_as=pin_file)
    print("\n")
    
    # merge features
    pin_file = os.path.basename(pos_file).split(".")[0].rsplit("_", 2)[0] + "_features"
    pin_file = '-'.join(filter(None, (pin_file, file_root))) + ".pin"
    pin_file = "{}/{}".format(output_dir, pin_file) if save_features == True else None
    df_merged = merge_feature_data(df_pos, df_neg, save_as=pin_file)

    print("\n\nDone!")


def format_to_pin(feature_file, data_type, save_as=None):
    '''Format feature files to pin.'''
    df = pd.read_csv(feature_file, sep="\t")
    print("Observations in {} data: {:,}".format(data_type, df.shape[0]))

    # add Label col of 0s
    df.insert(1, "Label", 0)
    
    # make Label value 1 for positive data and -1 for negative data
    label_val = 1 if data_type == "positive" else -1
    df["Label"] = df["Label"].replace(0, label_val)
    print("Data - rows: {:,} | columns: {:,}".format(df.shape[0], df.shape[1]))

    # remove non-feature columns
    to_remove = ["motif_name", "query_name"]
    df = df.drop(to_remove, axis=1)

    # create use motif_query column as SpecID
    df = df.rename(columns={"motif_query": "SpecID"})
    df["Peptide"] = df.index
    df["Proteins"] = df.index
    df["scannr"] = df.index
    print("Formatted {} data - rows: {:,} | columns: {:,}".format(data_type, df.shape[0], df.shape[1]))

    if save_as is not None:
        # save features to pin file
        df.to_csv(save_as, sep="\t", header=True, index=False)
        print("Saved formatted features to: {}".format(save_as))
    
    return df

def merge_feature_data(df1, df2, save_as):
    '''Merge pin formated feature files. (File must have the same column names.)'''
    # merge datasets
    df = pd.concat([df1, df2], axis=0)
    df = df.reset_index(drop=True)

    # reindex & copy to Peptide, Proteins, and scannr columns
    df["Peptide"] = df.index
    df["Proteins"] = df.index
    df["scannr"] = df.index
    print("Formatted merged data - rows: {:,} | columns: {:,}".format(df.shape[0], df.shape[1]))
    if (df1.shape[0] + df2.shape[0]) != df.shape[0]:
        raise Exception("Merged data rows not equal to input.\nmerged: {}, df1: {}, df2: {}".format(df.shape[0], df1.shape[0], df2.shape[0]))
    # remove rows with NaNs
    print("Total NaN values: {:,}".format(df.isna().sum().sum()))
    df = df.dropna()

    if save_as is not None:
        # save merged features to pin file
        df.to_csv(save_as, sep="\t", header=True, index=False)
        print("Saved merged features to: {}".format(save_as))

    return df


if __name__ == "__main__":
    
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pos_features", required=True, help="Positive feature data.")
    parser.add_argument("-n", "--neg_features", required=True, help="Negative feature data.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory.")
    parser.add_argument("-r", "--file_root", required=False, default=None, 
                        help="File root. Suffix will be appended to the END of the pin file name. (after features)")
    parser.add_argument("-s", "--save_features", required=False, default=True, help="Save formatted pin to output file.")
    args = parser.parse_args()
    
    # run script
    main(args)