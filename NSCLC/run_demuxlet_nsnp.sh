#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "ERROR: please provide num SNPs"
    exit 1
fi

nsnp=$1

source ~/miniforge3/bin/activate popscle
bam="40k_all.bam"

{ time ./run_demuxlet.sh $bam ${nsnp}.vcf.gz demuxlet_donor_all_${nsnp} 2> demuxlet_donor_all_${nsnp}.stderr; } 2> demuxlet_donor_all_${nsnp}.time

