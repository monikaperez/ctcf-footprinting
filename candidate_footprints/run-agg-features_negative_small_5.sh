#!/usr/bin/env bash
# Standard output and error log
#SBATCH --output=/mmfs1/gscratch/stergachislab/mwperez/logs/run-agg-features-negative_100bp-small_5_L_%j.out
#SBATCH -t 5-0
#SBATCH -c 5
#SBATCH --mem=80G

pwd; hostname; date
# once there is an error in the bash script the script ends with an error status
set -euo pipefail

# Aggregate features of 5% negative dataset

# Get fiberseq reads around CTCF motifs (using 5% of the data)
# use 100 bp flanking motif start

#echo "running agg_features.py for negative small 5% data set"
#python agg_features.py -i CTCF_negative_m6a_fiberseq_100bp_small_5.txt -o ../feature_data/

echo "running agg_features.py for L negative small 5% data set"
python agg_features.py -i CTCF_m6a_fiberseq_L_100bp_small_5_negative.txt -o ../feature_data/

echo -e "\n"
echo "DONE!"
