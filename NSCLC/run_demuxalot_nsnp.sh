#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "ERROR: please provide num SNPs"
    exit 1
fi

nsnp=$1

source /home/nkschaefer/miniforge3/etc/profile.d/conda.sh
conda activate demuxalot

bam="40k_all.bam"
donor="donor_all"

{ time ./run_demuxalot.py $bam donor_all.bcs ${nsnp}.vcf.gz 1> demuxalot_donor_all_${nsnp}.assignments 2> demuxalot_donor_all_${nsnp}.stderr; } 2> demuxalot_donor_all_${nsnp}.time

