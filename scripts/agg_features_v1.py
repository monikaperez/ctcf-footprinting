'''Aggregate features for CTCF model.'''

import os
import argparse

import pandas as pd
import numpy as np
from collections import defaultdict
from collections import Counter
import itertools
from functools import reduce

def main(args):

    # ensure input file exists
    input_file = os.path.abspath(args.input_file)
    output_dir = os.path.abspath(args.output_dir)
    assert os.path.exists(input_file), f"ERROR: input file not found: {input_file}"
    assert os.path.exists(output_dir), f"ERROR: output directory not found: {output_dir}"

    # read in args
    file_root = args.file_root
    agg_rle = args.rle

    # read in data
    df = pd.read_csv(input_file, sep="\t")
    print("ft center data - rows: {:,}, columns: {:,}".format(df.shape[0], df.shape[1]))

    # filter for only m6a & msp rows
    df = df[df["centered_position_type"].isin(["m6a", "msp"])]
    print("MSP & m6a instances - rows: {:,}, columns: {:,}".format(df.shape[0], df.shape[1]))

    # add columns of unique motif names & motif_query names
    make_IDs(df, col_name="motif_name")
    print_df_info(df)

    # remove rows with Ns in sequence
    df = clean_sequences(df)
    print("Cleaned rows: {:,}".format(df.shape[0]))
    
    # filter for regions within an MSP & add col with MSP size
    df = filt_msps(df)
    # remove MSP rows
    df = df[df["centered_position_type"] == "m6a"]
    print("Total m6a observations: " + "{:,}".format(df.shape[0]))
    # filter for m6a's within 40 bp range
    m6a_range_mask = (df["centered_start"] >= -40) & (df["centered_end"] <= 75)
    print("Total m6a's within 40 bp flank of motif: " + "{:,}".format(m6a_range_mask.sum()))
    df = df[m6a_range_mask]

    # group by motif/query name
    grouping_cols = ["motif_name", "query_name"]
    df_grouped = df.groupby(grouping_cols)
    # get group names (keys)
    group_names = list(df_grouped.groups.keys())
    print("Unique motif-sequence groups: " + "{:,}".format(len(group_names)))

    # extract features by motif/query group
    print_df_info(df)

    if agg_rle == False:
        # get all possible k-mers
        all_kmers = []
        for seq in df["subset_sequence"]:
            all_kmers.extend(get_kmers(get_motif_seq(seq), 3))
        # sorted unique k-mers containing AT
        all_kmers = [*set(all_kmers)]
        all_kmers = sorted([kmer for kmer in all_kmers if ("T" in kmer) or ("A" in kmer)])
        # make k-mer columns
        kmer_cols = list(itertools.chain(*[[kmer+"_count", kmer+"_m6a_prop"] for kmer in all_kmers]))

        # aggregate features
        print("\nAggregating features!")
        res = df.groupby(grouping_cols).apply(lambda x: agg_features(x, kmer_cols)).reset_index()
        # change col dtypes for memory (~50%)
        res[res.columns.tolist()[2:]] = res[res.columns.tolist()[2:]].apply(pd.to_numeric, downcast="float")

    else:
        print("Using rle features!")
        rle_ranges = [[0], [1], [2], [3], [4,5], [6,7], [8,10], [11,15], [16,20], [21,35], [5,35]]
        print("Using rle ranges:\n{}".format(["rle_" + "_".join(str(x) for x in rle_range) for rle_range in rle_ranges]))
        res = df.groupby(grouping_cols).apply(lambda x: agg_features_rle(x, rle_ranges)).reset_index()
        # change col dtypes for memory (~50%)
        res[res.columns.tolist()[2:]] = res[res.columns.tolist()[2:]].apply(pd.to_numeric, downcast="float")

    print("Total rows: " + "{:,}".format(res.shape[0]))
    print("Total columns: " + "{:,}".format(res.shape[1]))
    # add uniuqe motif/query name column
    res.insert(loc=0, column="motif_query", value=res["motif_name"].astype(str)+"/"+res["query_name"].astype(str))

    if agg_rle == False:
        print("Saving as file w/ rle.")
        output_file = os.path.join(output_dir, os.path.basename(input_file).replace(".txt", 
                                                                                    '-'.join(filter(None, ("_features", file_root))) + ".txt"))
    else:
        print("Saving as rle file.")
        output_file = os.path.join(output_dir, os.path.basename(input_file).replace(".txt", 
                                                                                    '-'.join(filter(None, ("_features-rle", file_root))) + ".txt"))
    print("Saving to output file: {}".format(output_file))
    res.to_csv(output_file, header=True, index=None, sep="\t",)

     
    print("Done!")


def clean_sequences(df):
    ''''Remove rows with an N character in sequence.'''
    print("Removing {:,}".format(df["subset_sequence"].str.contains("N").sum()))
    return df[~df["subset_sequence"].str.contains("N")]


def make_IDs(df, col_name="motif_name"):
    '''Makes a column of columns joined by a string.'''
    # add column to df of unique motif names
    df.insert(loc=0, column=col_name, value=df["chrom"]+"_"+df["centering_position"].astype(str)+"_"+df["strand"].astype(str))

def weird_division(n, d):
    '''Returns 0 if trying to divide by 0'''
    return n / d if d else 0

def filt_msps(df, motif_len=35):
    '''Filters for motif/query instances with a motif within 
    a MSP and adds MSP length to each motif/query group.'''
    
    # position of MSPs containing a motif
    msp_mask = (df["centered_position_type"] == "msp") & (df["centered_start"] <= 0) & (df["centered_end"] >= motif_len)
    print("MSP's with a motif: " + "{:,}".format(msp_mask.sum()))
    msp_groups = [(row["motif_name"], row["query_name"]) for idx, row in df[msp_mask].iterrows()]
    
    # filter for rows with motifs within an MSP (gets both MSP's & m6a's)
    df = df[df[["motif_name", "query_name"]].apply(tuple, 1).isin(msp_groups)]
    print("Total observations from motifs within an MSP: " + "{:,}".format(df.shape[0]))
    
    # add MSP size corresponding to each group
    # df with msp sizes
    df_msp = df.loc[msp_mask, ["motif_name", "query_name", "centered_end", "centered_start"]]
    df_msp["msp_size"] = df_msp["centered_end"] - df_msp["centered_start"]

    # match msp back to it's motif & fiber
    df = df.merge(df_msp[["motif_name", "query_name", "msp_size"]], on=["motif_name", "query_name"])
    print("Merged df shape w/ msp_size: " + "{:,}".format(df.shape[0]))
    return df


def print_df_info(df):
    '''Print df content and shape info.'''
    print("-----")
    print("Columns: {:,}".format(df.shape[1]))
    print("Total rows: " + "{:,}".format(df.shape[0]))
    print("Total unique motifs: {:,}".format(df.motif_name.nunique()))
    print("Total unique query_names: " + "{:,}".format(df.query_name.nunique()))
    print("MSP's: " + "{:,}".format(df[df["centered_position_type"] == "msp"].shape[0]))
    print("m6a's: " + "{:,}".format(df[df["centered_position_type"] == "m6a"].shape[0]))
    print("-----")


def get_motif_seq(subset_sequence):
    '''Get motif sequence from subset_sequence.'''
    center_idx = 100
    motif_len = 35
    motif_seq = subset_sequence[center_idx:(center_idx+motif_len)]
    return motif_seq


def get_kmers(seq, k):
    '''Decompose the input sequence (str) into k-mers (str).'''
    num_kmers = len(seq) - k + 1
    seqs = [seq[i:i+k] for i in range(num_kmers)]
    return seqs


def agg_features(x, feature_cols, rle_ranges):
    '''Collects features. Returns a df with feature info per unique group.'''
    feature_cols = ["msp_size",
                    "left_m6a_count", "right_m6a_count", "motif_m6a_count", # m6a count
                    "left_AT_count", "right_AT_count", "motif_AT_count", # AT count
                    "left_AT_prop", "right_AT_prop", "motif_AT_prop", # proportion of bases that are AT
                    "left_m6a_prop", "right_m6a_prop", "motif_m6a_prop", # proportion of ATs that are methylated
                    "total_m6a_prop", # proportion of ATs that are methylated in subseq_sequence
                   ] + feature_cols
    # add rle_range cols
    feature_cols = feature_cols + ["rle_" + "_".join(str(x) for x in rle_range) for rle_range in rle_ranges]
    
    # dict to hold info
    d = dict((col, 0) for col in feature_cols)
    
    # msp size
    d[feature_cols[0]] = x["msp_size"].unique()[0]
    
    # m6a count flanking left, right, motif
    d[feature_cols[1]] = (x["centered_end"] < 0).sum()
    d[feature_cols[2]] = (x["centered_start"] >= 35).sum()
    d[feature_cols[3]] = ((x["centered_start"] >= 0) & (x["centered_end"] <= 35)).sum()
    
    # sequences
    center = 100
    motif_len = 35
    subseq = x["subset_sequence"].unique()[0]
    subseq_motif = subseq[center:(center+motif_len)]
    # length of flank in bp
    flank_len = 40
    subseq_l = subseq[center-flank_len:center]
    subseq_r = subseq[center+motif_len:center+motif_len+flank_len]
    
    # AT count
    d[feature_cols[4]] = (subseq_l.count("A") + subseq_l.count("T"))
    d[feature_cols[5]] = (subseq_r.count("A") + subseq_r.count("T"))
    d[feature_cols[6]] = (subseq_motif.count("A") + subseq_motif.count("T"))
    
    # proportion of bases that are AT
    d[feature_cols[7]] = weird_division((subseq_l.count("A") + subseq_l.count("T")), flank_len)
    d[feature_cols[8]] = weird_division((subseq_r.count("A") + subseq_r.count("T")), flank_len)
    d[feature_cols[9]] = weird_division((subseq_motif.count("A") + subseq_motif.count("T")), motif_len)
    
    # proportion of methylated ATs (m6a_count/AT_count) return 0 if no ATs
    d[feature_cols[10]] = weird_division(d[feature_cols[1]], d[feature_cols[4]])
    d[feature_cols[11]] = weird_division(d[feature_cols[2]], d[feature_cols[5]])
    d[feature_cols[12]] = weird_division(d[feature_cols[3]], d[feature_cols[6]])
    # proportion of methylated ATs in total subset_sequence
    d[feature_cols[13]] = weird_division((x["centered_position_type"] == "m6a").sum(), (subseq.count("A") + subseq.count("T")))
    
    #----- k-mer info -----#
    # m6a instances in motif
    m6a_bool = x[(x["centered_start"] >= 0) & (x["centered_end"] <= 35)]["centered_start"].values
    m6a_bool = [1 if (i in m6a_bool) else 0 for i in range(0, 35)]
    
    # decompose motif seq into k-mers
    motif_kmers = get_kmers(subseq_motif, 3)
    m6a_kmers = get_kmers(m6a_bool, 3)
    
    # make bool of ATs in motif k-mers
    motifis_AT = [("A" in kmer or "T" in kmer) for kmer in motif_kmers]
    
    # drop motif instances without an AT
    motif_kmers = list(np.array(motif_kmers)[motifis_AT])
    m6a_kmers = list(np.array(m6a_kmers)[motifis_AT])
    
    # make kmer_dict
    kmer_info = defaultdict(lambda: defaultdict(int))
    # collect count & m6a count per k-mer
    for kmer, m6as in zip(motif_kmers, m6a_kmers):
        kmer_info[kmer]["count"] += 1
        kmer_info[kmer]["m6a"] += sum(m6as)
    
    # get kmer_m6a proportion
    for kmer, v in kmer_info.items():
        d[kmer+"_count"] = v["count"]
        d[kmer+"_m6a_prop"] = weird_division(v["m6a"], ((Counter(kmer)["T"] + Counter(kmer)["A"]) * v["count"]))


    #----- rle info -----#

    # AT mask within motif sequence
    motifis_AT = [(base == "A" or base == "T") for base in subseq_motif]
    # take motif m6a & subset by positions that are actually AT
    m6as = np.array(m6a_bool)[np.array(motifis_AT)]
    
    # get rle counts in ranges
    if 0 in m6as:
        # get rle
        rle_res = rle(m6as)
        # add run_lengths of 0 (adjacent m6a's)
        d["rle_0"] = reduce(lambda sum, j: sum + (j-1 if j > 1 else 0), rle_res[0][rle_res[2] == 1], 0)
        # get run_lengths when value = 0
        rle_res = rle_res[0][rle_res[2] == 0]
    
        # count number of instances within ranges
        for rle_range in rle_ranges[1:]:
            col_name = "rle_" + "_".join(str(x) for x in rle_range)
            d[col_name] = count_range_in_list(rle_res, *rle_range)
    else:
        for rle_range in rle_ranges:
            col_name = "rle_" + "_".join(str(x) for x in rle_range)
            d[col_name] = 0
    
    return pd.Series(d, index=list(d.keys()))

def rle(inarray):
    ''' run length encoding. Partial credit to R rle function. 
        Multi datatype arrays catered for including non Numpy
        returns: tuple (runlengths, startpositions, values) 
        input: array
        output: tuple of (run_lengths, start_positions, values)'''
    ia = np.asarray(inarray)                # force numpy
    n = len(ia)
    if n == 0: 
        return (None, None, None)
    else:
        y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)   # must include last element position
        z = np.diff(np.append(-1, i))       # run lengths
        p = np.cumsum(np.append(0, z))[:-1] # positions
        return(z, p, ia[i])
    
def count_range_in_list(li, min, max=None):
    '''Count instances of values within range in list.'''
    max = min if not max else max
    ctr = 0
    for x in li:
        if min <= x <= max:
            ctr += 1
    return ctr

def agg_features_rle(x, rle_ranges):
    '''Collects features with rle! Returns a df with feature info per unique group.'''
    feature_cols = ["msp_size",
                    "left_m6a_count", "right_m6a_count", "motif_m6a_count", # m6a count
                    "left_AT_count", "right_AT_count", "motif_AT_count", # AT count
                    "left_AT_prop", "right_AT_prop", "motif_AT_prop", # proportion of bases that are AT
                    "left_m6a_prop", "right_m6a_prop", "motif_m6a_prop", # proportion of ATs that are methylated
                    "total_m6a_prop", # proportion of ATs that are methylated in subseq_sequence
                   ]
    # add rle_range cols
    feature_cols = feature_cols + ["rle_" + "_".join(str(x) for x in rle_range) for rle_range in rle_ranges]
    
    # dict to hold info
    d = dict((col, 0) for col in feature_cols)
    
    # msp size
    d[feature_cols[0]] = x["msp_size"].unique()[0]
    
    # m6a count flanking left, right, motif
    d[feature_cols[1]] = (x["centered_end"] < 0).sum()
    d[feature_cols[2]] = (x["centered_start"] >= 35).sum()
    d[feature_cols[3]] = ((x["centered_start"] >= 0) & (x["centered_end"] <= 35)).sum()
    
    # sequences
    center = 100
    motif_len = 35
    subseq = x["subset_sequence"].unique()[0]
    subseq_motif = subseq[center:(center+motif_len)]
    # length of flank in bp
    flank_len = 40
    # flank left & flank_right
    subseq_l = subseq[center-flank_len:center]
    subseq_r = subseq[center+motif_len:center+motif_len+flank_len]
    
    # AT count
    d[feature_cols[4]] = (subseq_l.count("A") + subseq_l.count("T"))
    d[feature_cols[5]] = (subseq_r.count("A") + subseq_r.count("T"))
    d[feature_cols[6]] = (subseq_motif.count("A") + subseq_motif.count("T"))
    
    # proportion of bases that are AT
    d[feature_cols[7]] = (subseq_l.count("A") + subseq_l.count("T"))/flank_len
    d[feature_cols[8]] = (subseq_r.count("A") + subseq_r.count("T"))/flank_len
    d[feature_cols[9]] = (subseq_motif.count("A") + subseq_motif.count("T"))/motif_len
    
    # proportion of methylated ATs (m6a_count/AT_count)
    d[feature_cols[10]] = d[feature_cols[1]]/d[feature_cols[4]]
    d[feature_cols[11]] = d[feature_cols[2]]/d[feature_cols[5]]
    d[feature_cols[12]] = d[feature_cols[3]]/d[feature_cols[6]]
    d[feature_cols[13]] = (x["centered_position_type"] == "m6a").sum()/(subseq.count("A") + subseq.count("T"))
    
     #----- rle info -----#

    # m6a instances in motif
    m6a_bool = x[(x["centered_start"] >= 0) & (x["centered_end"] <= 35)]["centered_start"].values
    m6a_bool = [1 if (i in m6a_bool) else 0 for i in range(0, 35)]
    # AT mask
    motifis_AT = [(base == "A" or base == "T") for base in subseq_motif]
    # take motif m6a & subset by positions that are actually AT
    m6as = np.array(m6a_bool)[np.array(motifis_AT)]
    
    # get rle counts in ranges
    if 0 in m6as:
        # get rle
        rle_res = rle(m6as)
        # add run_lengths of 0 (adjacent m6a's)
        d["rle_0"] = reduce(lambda sum, j: sum + (j-1 if j > 1 else 0), rle_res[0][rle_res[2] == 1], 0)
        # get run_lengths when value = 0
        rle_res = rle_res[0][rle_res[2] == 0]
    
        # count number of instances within ranges
        for rle_range in rle_ranges[1:]:
            col_name = "rle_" + "_".join(str(x) for x in rle_range)
            d[col_name] = count_range_in_list(rle_res, *rle_range)
    else:
        for rle_range in rle_ranges:
            col_name = "rle_" + "_".join(str(x) for x in rle_range)
            d[col_name] = 0
    
    return pd.Series(d, index=list(d.keys()))



if __name__ == "__main__":
    
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True, help="ft center output file.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory to save feature file to.")
    parser.add_argument("-r", "--rle", required=False, default=False, help="Get rle features without k-mers.")
    parser.add_argument("-fr", "--file_root", required=False, default=None, 
                    help="File root. Suffix will be appended to the END of the pin file name. (after features)")
    args = parser.parse_args()
    
    # run script
    main(args)