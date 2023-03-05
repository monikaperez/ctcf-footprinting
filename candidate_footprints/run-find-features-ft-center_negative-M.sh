#!/usr/bin/env bash
# Standard output and error log
#SBATCH --output=/mmfs1/gscratch/stergachislab/mwperez/logs/run-find-features-ft-center-negative_100bp-M_%j.out
#SBATCH -t 5-0

pwd; hostname; date
# once there is an error in the bash script the script ends with an error status
set -euo pipefail

# Generating NEGATIVE dat aset

# Get fiberseq reads around CTCF motifs
# use 100 bp flanking motif start
# filter for m6A type
echo "running negative M motifs"
ft center -d 100 ../data/GM12878.aligned.bam <(cut -f 1,2,3,6 CTCF_negative_candidate_footprints_M.bed) > CTCF_negative_m6a_fiberseq_M_100bp.txt

echo -e "\n"
echo "DONE!"
