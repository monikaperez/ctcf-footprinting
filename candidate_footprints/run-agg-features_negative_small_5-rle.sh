#!/usr/bin/env bash
# Standard output and error log
#SBATCH --output=/mmfs1/gscratch/stergachislab/mwperez/logs/run-agg-features-L_100bp_positive-rle_%j.out
#SBATCH -t 5-0
#SBATCH -c 15
#SBATCH --mem=80G

pwd; hostname; date
# once there is an error in the bash script the script ends with an error status
set -euo pipefail

# Aggregate rle features of 5% negative dataset
#fiberseq_data="CTCF_m6a_fiberseq_L_100bp_small_5_negative.txt"
#echo "running agg_features.py for L negative small 5% data set (RLE)"
#echo "Using $fiberseq_data"
#python agg_features.py -i "$fiberseq_data" -o ../feature_data/ --rle True
#echo -e "\n"
#echo "DONE!"

# Aggregate rle features of positive dataset
fiberseq_data="CTCF_m6a_fiberseq_L_100bp_positive.txt"
echo "running agg_features.py for L POSITIVE data set (RLE)"
echo "Using $fiberseq_data"
python agg_features.py -i "$fiberseq_data" -o ../feature_data/ --rle True

echo -e "\n"
echo "DONE!"