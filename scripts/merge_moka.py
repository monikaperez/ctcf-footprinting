'''Merge mokapot FDR value with fiberseq data. Also saves fiberseq data filtered & formatted. '''

import os
import argparse
import pandas as pd

def main(args):

    # ensure input files exists
    input_file = os.path.abspath(args.input_file)
    pos_file = os.path.abspath(args.pos_file)
    neg_file = os.path.abspath(args.neg_file)
    output_dir = os.path.abspath(args.output_dir)
    assert os.path.exists(input_file), f"ERROR: file not found: {input_file}"
    assert os.path.exists(pos_file), f"ERROR: file not found: {pos_file}"
    assert os.path.exists(neg_file), f"ERROR: file not found: {neg_file}"
    assert os.path.exists(output_dir), f"ERROR: output directory not found: {output_dir}"

    print("Mokapot res file: {}".format(os.path.basename(input_file)))
    print("Positive data file: {}".format(os.path.basename(pos_file)))
    print("Negative data file: {}".format(os.path.basename(neg_file)))
    print("Saving output to: {}".format(output_dir))
    print("\n")

    # read in args
    file_root = args.file_root
    #save_features = args.save_features

    # read mokapot res
    res = pd.read_csv(input_file, sep="\t")
    print("mokapot data - rows: {:,} | cols: {:,}".format(res.shape[0], res.shape[1]))
    # rename columns to match fiberseq data
    res = res.rename(columns={"SpecID": "motif_query", "mokapot q-value": "FDR"})

    # ----- merge positive data -----
    # read fiberseq file
    print("Merging positive data.")
    d_p = pd.read_csv(pos_file, sep="\t")
    # merge w/ mokapot data
    print("data - rows: {:,} | cols: {:,}".format(d_p.shape[0], d_p.shape[1]))
    d_p = pd.merge(d_p, res[["motif_query", "Label", "FDR"]], on="motif_query", how="inner")
    print("merged - rows: {:,} | cols: {:,}".format(d_p.shape[0], d_p.shape[1]))
    # change Label from True to 1
    d_p.Label = d_p.Label.replace({True: 1})
    # save to file
    output_file = os.path.basename(input_file).replace("psms", "m6a_fiberseq-positive")
    output_file = os.path.join(output_dir, output_file)
    d_p.to_csv(output_file, sep="\t", header=True, index=False)
    print("Saving positive data to: {}".format(os.path.basename(output_file)))
    
    # ----- merge negative data -----
    print("Formatting negative data.")
    d_n = pd.read_csv(neg_file, sep="\t")
    print("data - rows: {:,} | cols: {:,}".format(d_n.shape[0], d_n.shape[1]))
    # add Label & fake mokapot q-value/FDR to NEG table
    d_n["Label"] = -1
    d_n["FDR"] = 1
    print("merged - rows: {:,} | cols: {:,}".format(d_n.shape[0], d_n.shape[1]))
    # save to file
    output_file = os.path.basename(input_file).replace("psms", "m6a_fiberseq-negative")
    output_file = os.path.join(output_dir, output_file)
    d_n.to_csv(output_file, sep="\t", header=True, index=False)
    print("Saving negative data to: {}".format(os.path.basename(output_file)))

    # ----- merge pos & neg data -----
    print("Merging positive and negative data.")
    print("Positive data - rows: {:,} | cols: {:,}".format(d_p.shape[0], d_p.shape[1]))
    print("Negative data - rows: {:,} | cols: {:,}".format(d_n.shape[0], d_n.shape[1]))
    assert d_p.columns.tolist() == d_n.columns.tolist(), f"ERROR: positive and negative data do not have the same columns."
    # merge negative features w/ positive features (append rows)
    d_m = pd.concat([d_p, d_n])
    print("merged data - rows: {:,} | cols: {:,}".format(d_m.shape[0], d_m.shape[1]))
    # save merged table
    output_file = os.path.basename(input_file).replace("psms", "m6a_fiberseq")
    output_file = os.path.join(output_dir, output_file)
    d_m.to_csv(output_file, sep="\t", header=True, index=False)
    print("Saving merged file to: {}".format(os.path.basename(output_file)))

if __name__ == "__main__":
    
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True, help="Output mokapot file.")
    parser.add_argument("-p", "--pos_file", required=True, help="Cleaned positive m6a fiberseq data.")
    parser.add_argument("-n", "--neg_file", required=True, help="Cleaned negative m6a fiberseq data.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory.")
    parser.add_argument("-r", "--file_root", required=False, default=None, 
                        help="File root. Suffix will be appended to the END of the pin file name. (after features)")
    #parser.add_argument("-s", "--save_features", required=False, default=True, help="Save formatted pin to output file.")
    args = parser.parse_args()
    
    # run script
    main(args)