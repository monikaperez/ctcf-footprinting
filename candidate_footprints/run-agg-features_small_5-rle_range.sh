#!/usr/bin/env bash
# Standard output and error log
#SBATCH --output=/mmfs1/gscratch/stergachislab/mwperez/logs/run-agg-features-L_100bp_small_5-rle_v2_negative_%j.out
#SBATCH -t 5-0
#SBATCH -c 15
#SBATCH --mem=80G

pwd; hostname; date
# once there is an error in the bash script the script ends with an error status
set -euo pipefail


# ----- POSITIVE -----

# Aggregate rle features of positive dataset
# specific rle ranges

#fiberseq_data="CTCF_m6a_fiberseq_L_100bp_positive.txt"
#echo "running agg_features.py for L POSITIVE data set (RLE)"
#echo "Using $fiberseq_data"
#
#python ../scripts/agg_features.py -i "$fiberseq_data" -o ../feature_data/ -fr v2
#
#echo -e "\n"
#echo "DONE!"


# Aggregate rle features of positive dataset
# getting ALL rle's from 0-35
#fiberseq_data="CTCF_m6a_fiberseq_L_100bp_positive.txt"
#echo "running agg_features.py for L POSITIVE data set (RLE)"
#echo "Getting ALL rle bins (0-35)."
#echo "Using $fiberseq_data"
#
#python ../scripts/agg_features.py -i "$fiberseq_data" -o ../feature_data/ -b all -fr v3
#
#echo -e "\n"
#echo "DONE!"


# Aggregate rle features of positive dataset
# getting ALL rle's from 0-35 & rle_max
fiberseq_data="CTCF_m6a_fiberseq_L_100bp_positive.txt"
echo "running agg_features.py for L POSITIVE data set (RLE)"
echo "Getting ALL rle bins (0-35)."
echo "Using $fiberseq_data"

python ../scripts/agg_features.py -i "$fiberseq_data" -o ../feature_data/ -b all -fr v4

echo -e "\n"
echo "DONE!"


# ----- NEGATIVE -----

# Aggregate rle features of 5% negative dataset
#fiberseq_data="CTCF_m6a_fiberseq_L_100bp_small_5_negative.txt"
#echo "running agg_features.py for L negative small 5% data set (RLE)"
#echo "Using $fiberseq_data"

#python ../scripts/agg_features.py -i "$fiberseq_data" -o ../feature_data/ -fr v2

#echo -e "\n"
#echo "DONE!"

# Aggregate rle features of positive dataset
# getting ALL rle's from 0-35
#fiberseq_data="CTCF_m6a_fiberseq_L_100bp_small_5_negative.txt"
#echo "running agg_features.py for L NEGATIVE data set (RLE)"
#echo "Getting ALL rle bins (0-35)."
#echo "Using $fiberseq_data"
#
#python ../scripts/agg_features.py -i "$fiberseq_data" -o ../feature_data/ -b all -fr v3
#
#echo -e "\n"
#echo "DONE!"


# Aggregate rle features of positive dataset
# getting ALL rle's from 0-35 and rle_max
#fiberseq_data="CTCF_m6a_fiberseq_L_100bp_small_5_negative.txt"
#echo "running agg_features.py for L NEGATIVE data set (RLE)"
#echo "Getting ALL rle bins (0-35)."
#echo "Using $fiberseq_data"

#python ../scripts/agg_features.py -i "$fiberseq_data" -o ../feature_data/ -b all -fr v4

#echo -e "\n"
#echo "DONE!"
