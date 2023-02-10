#!/usr/bin/env bash
# Standard output and error log
#SBATCH --output=/mmfs1/gscratch/stergachislab/mwperez/logs/run-find-features-ft-center-100bp_%j.out

pwd; hostname; date
# once there is an error in the bash script the script ends with an error status
set -euo pipefail

# Get fiberseq reads around CTCF motifs
# use 100 bp flanking motif start
# filter for m6A type (all types: m6a, 5mC, nuc, msp)
# gets sequence x bps around motif start
echo "running M motifs"
ft center -d 100 ../data/GM12878.aligned.bam <(cut -f 1,2,3,6 CTCF_candidate_footprints_M.bed) > CTCF_m6a_fiberseq_M_100bp.txt
echo "running L motifs"
ft center -d 100 ../data/GM12878.aligned.bam <(cut -f 1,2,3,6 CTCF_candidate_footprints_L.bed) > CTCF_m6a_fiberseq_L_100bp.txt
echo "running XL motifs"
ft center -d 100 ../data/GM12878.aligned.bam <(cut -f 1,2,3,6 CTCF_candidate_footprints_XL.bed) > CTCF_m6a_fiberseq_XL_100bp.txt
