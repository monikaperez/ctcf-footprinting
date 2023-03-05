#!/usr/bin/env bash
# Standard output and error log
#SBATCH --output=/mmfs1/gscratch/stergachislab/mwperez/logs/run-find-features-ft-center-negative_100bp-small_5_L_%j.out
#SBATCH -t 5-0
#SBATCH -c 5

pwd; hostname; date
# once there is an error in the bash script the script ends with an error status
set -euo pipefail

# Generating NEGATIVE dat aset

# Get fiberseq reads around CTCF motifs (using 1% of the data)
# use 100 bp flanking motif start
# filter for m6A type

#echo "running negative M motifs (small)"
#ft center -d 100 ../data/GM12878_small.aligned.bam <(cut -f 1,2,3,6 CTCF_negative_candidate_footprints_M.bed) > CTCF_negative_m6a_fiberseq_M_100bp_small.txt
#echo "running negative L motifs (small)"
#ft center -d 100 ../data/GM12878_small.aligned.bam <(cut -f 1,2,3,6 CTCF_negative_candidate_footprints_L.bed) > CTCF_negative_m6a_fiberseq_L_100bp_small.txt
#echo "running negative XL motifs (small)"
#ft center -d 100 ../data/GM12878_small.aligned.bam <(cut -f 1,2,3,6 CTCF_negative_candidate_footprints_XL.bed) > CTCF_negative_m6a_fiberseq_XL_100bp_small.txt

#echo "running merged negative motifs (small)"
#ft center -d 100 ../data/GM12878_small_5.aligned.bam <(cut -f 1,2,3,6 CTCF_negative_candidate_footprints.bed) > CTCF_negative_m6a_fiberseq_100bp_small_5.txt

echo "running L negative motifs (small 5%)"
ft center -d 100 ../data/GM12878_small_5.aligned.bam <(cut -f 1,2,3,6 CTCF_candidate_footprints_L_negative.bed) > CTCF_m6a_fiberseq_L_100bp_small_5.txt

echo -e "\n"
echo "DONE!"
