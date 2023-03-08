'''Aggregate features for CTCF model.'''

import os
import argparse

import pandas as pd
import numpy as np
from collections import defaultdict
from collections import Counter
import itertools
from functools import reduce
import logging

def main(args):

    # ensure input file exists
    input_file = os.path.abspath(args.input_file)
    output_dir = os.path.abspath(args.output_dir)
    assert os.path.exists(input_file), f"ERROR: input file not found: {input_file}"
    assert os.path.exists(output_dir), f"ERROR: output directory not found: {output_dir}"

    # read in args
    file_root = args.file_root
    save_cleaned = args.save_cleaned
    #rle_ranges = args.bins

    # read in data
    df = pd.read_csv(input_file, sep="\t")
    print("ft center data - rows: {:,}, columns: {:,}".format(df.shape[0], df.shape[1]))

    # filter for only m6a & msp rows
    df = df[df["centered_position_type"].isin(["m6a", "msp"])]
    print("MSP & m6a instances - rows: {:,}, columns: {:,}".format(df.shape[0], df.shape[1]))

    # remove rows with Ns in sequence & clean chroms
    df = clean_sequences(df)
    df = clean_chroms(df)
    print("Cleaned rows: {:,}".format(df.shape[0]))

    # add columns of unique motif names & motif_query names
    df.insert(loc=0, column="motif_name", value=df["chrom"]+"_"+df["centering_position"].astype(str)+"_"+df["strand"].astype(str))
    df.insert(loc=0, column="motif_query", value=df["motif_name"].astype(str)+"/"+df["query_name"].astype(str))
    
    # filter for regions within an MSP, MSP size per motif_query, remove MSP size
    df = filt_msps(df)
    df = df[df["centered_position_type"] == "m6a"]
    print("Total m6a observations: " + "{:,}".format(df.shape[0]))

    # save filtered m6a fiberseq data (within 100bp)
    if save_cleaned:
        output_file = input_file.replace(".txt", "-cleaned-100bp.txt")
        print("Saving to cleaned m6a fiberseq data to: {}".format(output_file))
        df.to_csv(output_file, header=True, index=None, sep="\t",)

    # filter for m6a's within 40 bp range
    m6a_range_mask = (df["centered_start"] >= -40) & (df["centered_end"] <= 75)
    print("Total m6a's within 40 bp flank of motif: " + "{:,}".format(m6a_range_mask.sum()))
    df = df[m6a_range_mask]

    # group by motif/query name AND centered_query_start & centered_query_end 
    # (same read aligning twice, multiple repeats in a row double counting m6A instances but NOT sequences)
    grouping_cols = ["motif_name", "query_name"]
    df_grouped = df.groupby(grouping_cols)
    # get group names (keys)
    group_names = list(df_grouped.groups.keys())
    print("Unique motif-fiber groups: " + "{:,}".format(len(group_names)))

    # extract features by motif/query group
    print("\nAggregating features!")
    res = df.groupby(grouping_cols).apply(lambda x: agg_features(x)).reset_index()
    print("Features: {}".format(res.columns.tolist()[2:]))
    print("Total rows: " + "{:,}".format(res.shape[0]))
    print("Total columns: " + "{:,}".format(res.shape[1]))

    print("Removing rows w/ any proportion > 1. (no idea why.)")
    res = res[~(res.loc[:,res.columns.str.contains("prop")] > 1).any(1)]
    print("Total rows: " + "{:,}".format(res.shape[0]))

    # add uniuqe motif/query name column
    res.insert(loc=0, column="motif_query", value=res[grouping_cols].apply(lambda row: "/".join(row.values.astype(str)), axis=1))

    print("Saving with k-mer and rle file.")
    output_file = os.path.join(output_dir, os.path.basename(input_file).replace(".txt", '-'.join(filter(None, ("_features", file_root))) + ".txt"))
    print("Saving to output file: {}".format(output_file))
    res.to_csv(output_file, header=True, index=None, sep="\t",)
     
    print("Done!")


def clean_sequences(df):
    ''''Remove rows with an N character in sequence.'''
    print("Removing sequences w/ N: {:,}".format(df["subset_sequence"].str.contains("N").sum()))
    return df[~df["subset_sequence"].str.contains("N")]


def clean_chroms(df):
    ''''Remove rows no in standard chromosomes (chr1-22, chrX).'''
    clean_chroms = (df["chrom"].isin([f"chr{x}" for x in list(range(1, 23))+["X"]]))
    print("Removing weird chromosomes: {:,}".format((df.shape[0] - clean_chroms.sum())))
    return df.loc[clean_chroms]


def revcomp(seq):
    ''''Get the reverse complement of sequence.'''
    COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(COMPLEMENT[base] for base in reversed(seq))


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
    dup_alignments = df_msp.shape[0]
    df_msp = df_msp.drop_duplicates(subset=["motif_name", "query_name"], keep=False)
    print("Dropped motif-query instances w/ multiple MSPs (multiple alignments). \
        {:,}".format(dup_alignments - df_msp.shape[0]))

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


# make all possible 3-mers
def make_kmer_dict(kmer_size, use_canonical=False):
    kmer_dict = {}
    kmers = []
    # cartesian product (set formed from 2+ sets w/ all ordered pairs)
    for kmer_list in itertools.product(["A", "C", "G", "T"], repeat=kmer_size):
        kmer = "".join(kmer_list)
        rc_kmer = revcomp(kmer)
        if rc_kmer < kmer and use_canonical:
            kmer = rc_kmer
        # skip no m6a possible kmers
        if kmer.count("A") + kmer.count("T") == 0:
            continue
        kmers.append(kmer)
        kmer_dict[f"{kmer}_count"] = 0
        kmer_dict[f"{kmer}_m6a_count"] = 0
    return kmer_dict, kmers


def kmer_features(seq, m6a_bool, kmer_size = 3, use_canonical = False):
    stop_index = len(seq) - kmer_size + 1
    kmer_counts, kmers = make_kmer_dict(kmer_size, use_canonical=use_canonical)
    kmer_feats = {}
    for i in range(stop_index):
        # get current index
        kmer = seq[i:i+kmer_size]
        m6a_count = m6a_bool[i:i+kmer_size].sum()
        # skip no m6a possible kmers
        if kmer.count("A") + kmer.count("T") == 0:
            continue
        # get canonical kmer
        rc_kmer = revcomp(kmer)
        if rc_kmer < kmer and use_canonical:
            kmer = rc_kmer
        # counts 
        kmer_counts[f"{kmer}_count"] += 1
        kmer_counts[f"{kmer}_m6a_count"] += m6a_count
        
    # get motif_m6a_prop
    for kmer in kmers:
        AT_count = (kmer.count("A") + kmer.count("T"))
        motif_m6a_prop = weird_division(kmer_counts[f"{kmer}_m6a_count"], 
                                        (AT_count*kmer_counts[f"{kmer}_count"]))
        #if motif_m6a_prop > 1:
        #    raise ValueError(f"Motif m6A proportion > 1. {kmer} value: {motif_m6a_prop}")
        kmer_feats[f"{kmer}_count"] = kmer_counts[f"{kmer}_count"]
        kmer_feats[f"{kmer}_m6a_prop"] = weird_division(kmer_counts[f"{kmer}_m6a_count"],
                                                        (AT_count * kmer_counts[f"{kmer}_count"]))
    return kmer_feats


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


def agg_features(x):
    '''Collect motif features per group.'''
    d = defaultdict()
    
    # msp size
    d["msp_size"] = x["msp_size"].unique()[0]
    
    # get sequences (left, right, motif)
    center = 100
    motif_len = 35
    flank = 40
    subseq = x["subset_sequence"].unique()[0]
    subseq_motif = subseq[center:(center+motif_len)]
    subseq_l = subseq[(center-flank):center]
    subseq_r = subseq[(center+motif_len):(center+motif_len+flank)]
    
    # AT count (left, right, motif)
    d["left_AT_count"] = (subseq_l.count("A") + subseq_l.count("T"))
    d["right_AT_count"] = (subseq_r.count("A") + subseq_r.count("T"))
    d["motif_AT_count"] = (subseq_motif.count("A") + subseq_motif.count("T"))
    
    # proportion of bases that are AT
    d["left_AT_prop"] = weird_division(d["left_AT_count"], flank)
    d["right_AT_prop"] = weird_division(d["right_AT_count"], flank)
    d["motif_AT_prop"] = weird_division(d["motif_AT_count"], motif_len)
    
    # m6a instances
    d["left_m6a_count"] = (x["centered_start"] < 0).sum()
    d["right_m6a_count"] = (x["centered_start"] >= 35).sum()
    d["motif_m6a_count"] = ((x["centered_start"] >= 0) & (x["centered_start"] < 35)).sum()
    
    # proportion of ATs that are methylated (m6a_count/AT_count) return 0 if no ATs
    d["left_m6a_prop"] = weird_division(d["left_m6a_count"], d["left_AT_count"])
    d["right_m6a_prop"] = weird_division(d["right_m6a_count"], d["right_AT_count"])
    d["motif_m6a_prop"] = weird_division(d["motif_m6a_count"], d["motif_AT_count"])
    
    # sanity check proportions
    for i in ["left", "right", "motif"]:
        if (d[f"{i}_AT_prop"] > 1) or (d[f"{i}_m6a_prop"] > 1):
            print("motif: ", x["motif_name"].unique()[0])
            print("query name: ", x["query_name"].unique()[0])
            
            print(f"{i} AT prop: ", d[f"{i}_AT_prop"])
            print(f"{i} m6a prop: ", d[f"{i}_m6a_prop"])
            print(f"{i} m6a count: ", d[f"{i}_m6a_count"])
            print(f"{i} AT count: ", d[f"{i}_AT_count"])
            
            print("centered position type(s): {}".format(x.centered_position_type.value_counts()))
            logging.error("AT or m6a proportion > 1.")
    
    # ----- kmers -----
    # m6a instances in motif
    m6a_bool = x[(x["centered_start"] >= 0) & (x["centered_end"] <= 35)]["centered_start"].values
    m6a_bool = [1 if (i in m6a_bool) else 0 for i in range(0, 35)]
    # convert m6a bool to numpy array
    m6a_bool = np.array(m6a_bool)
    # get counts & m6a prop per k-mer
    kmer_counts = kmer_features(subseq_motif, m6a_bool, kmer_size=3, use_canonical=True)
    
    # ----- rle -----
    # make np array for sequence
    np_subseq_motif = np.frombuffer(bytes(subseq_motif, "utf-8"), dtype="S1")
    # make AT mask
    AT_bases = (np_subseq_motif == b"A" ) | (np_subseq_motif == b"T")
    # get run lengths using the AT base space
    run_lengths, run_starts, run_values = rle(m6a_bool[AT_bases])
    # subset the run lengths that encode zeros (non-m6a):
    non_m6a_run_lengths = run_lengths[run_values == 0]

    d["rle_max"] = max(non_m6a_run_lengths) if non_m6a_run_lengths.size else 0
    
    d = d | kmer_counts
    
    return pd.Series(d, index=list(d.keys()))


if __name__ == "__main__":
    
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True, help="ft center output file.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory to save feature file to.")
    parser.add_argument("-fr", "--file_root", required=False, default=None, 
                        help="File root. Suffix will be appended to the END of the pin file name. (after features)")
    #parser.add_argument("-b", "--bins", required=False, default="all", help="Set bins for rle. Deafult is all (bins 0-35)")
    parser.add_argument("-c", "--save_cleaned", required=False, default=True, help="Save cleaned fiber file.")    
    args = parser.parse_args()
    
    # run script
    main(args)