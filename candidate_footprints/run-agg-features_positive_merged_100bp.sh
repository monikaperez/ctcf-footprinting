#!/usr/bin/env bash
# Standard output and error log
#SBATCH --output=/mmfs1/gscratch/stergachislab/mwperez/logs/run-agg-features-positive-merged-100bp_%j.out
#SBATCH -t 5-0
#SBATCH -c 5
#SBATCH --mem=80G

pwd; hostname; date
# once there is an error in the bash script the script ends with an error status
set -euo pipefail

# Aggregate features of merged positive data

echo "Running agg_features.py over merged 100bp positive data set (43,760,948 rows)"
echo "Using file: CTCF_m6a_fiberseq_merged_100bp_positive.txt"
echo "Saving output to: ../feature_data/"
python agg_features.py -i CTCF_m6a_fiberseq_merged_100bp_positive.txt -o ../feature_data/

echo -e "\n"
echo "DONE!"
