#!/usr/bin/env bash
# Standard output and error log
#SBATCH --output=/mmfs1/gscratch/stergachislab/mwperez/logs/run-find-features-intersect-ft-center-negative-100bp_%j.out
#SBATCH -t 1-0

pwd; hostname; date
# once there is an error in the bash script the script ends with an error status
set -euo pipefail

# Generating NEGATIVE dat aset
# Get motifs from Merged_CTCF_motifs.bed that do NOT intersect with peaks in CTCF-ChIP-peaks-and-DNase-hotspot.bed
#bedtools intersect -v -a ../fimo/Merged_CTCF_motifs.bed -b ../data/CTCF-ChIP-peaks-and-DNase-hotspot.bed > CTCF_negative_candidate_footprints.bed

# separate by motif type (optional)
#grep CTCF_M CTCF_negative_candidate_footprints.bed > CTCF_negative_candidate_footprints_M.bed
#grep CTCF_L CTCF_negative_candidate_footprints.bed > CTCF_negative_candidate_footprints_L.bed
#grep CTCF_XL CTCF_negative_candidate_footprints.bed > CTCF_negative_candidate_footprints_XL.bed

# Get fiberseq reads around CTCF motifs
# use 100 bp flanking motif start
# filter for m6A type
echo "running negative M motifs"
ft center -d 100 ../data/GM12878.aligned.bam <(cut -f 1,2,3,6 CTCF_negative_candidate_footprints_M.bed) > CTCF_negative_m6a_fiberseq_M_100bp.txt
echo "running negative L motifs"
ft center -d 100 ../data/GM12878.aligned.bam <(cut -f 1,2,3,6 CTCF_negative_candidate_footprints_L.bed) > CTCF_negative_m6a_fiberseq_L_100bp.txt
echo "running negative XL motifs"
ft center -d 100 ../data/GM12878.aligned.bam <(cut -f 1,2,3,6 CTCF_negative_candidate_footprints_XL.bed) > CTCF_negative_m6a_fiberseq_XL_100bp.txt

echo "DONE!"
