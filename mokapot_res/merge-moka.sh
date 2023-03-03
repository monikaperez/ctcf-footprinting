#!/usr/bin/env bash
# Standard output and error log
#SBATCH --output=/mmfs1/gscratch/stergachislab/mwperez/logs/merge_moka_CTCF_100bp_small_5_L_v2_%j.out
#SBATCH -t 5-0
#SBATCH -c 10

pwd; hostname; date
# once there is an error in the bash script the script ends with an error status
set -euo pipefail
# get execution time
start=`date +%s`

#----- CTCF L v2 -----

# Merge mokapot res with fiberseq data
input_file="CTCF_100bp_L_5_v2.mokapot.psms.txt"
pos_file="CTCF_m6a_fiberseq_L_100bp_positive-cleaned_100bp.txt"
neg_file="CTCF_m6a_fiberseq_L_100bp_small_5_negative-cleaned_100bp.txt"
file_root="100bp_L_5_v2"

echo "Merging mokapot results"
echo "Mokapot file: $input_file"
echo "Positive file: $pos_file"
echo "Negative file: $neg_file"
echo "File root: $file_root"
echo -e "\n"
python ../scripts/merge_moka.py -i $input_file -p $pos_file -n $neg_file -o . -r $file_root

echo -e "\n"
echo "DONE!"



end=`date +%s`
runtime=$((end-start))
echo "$runtime"