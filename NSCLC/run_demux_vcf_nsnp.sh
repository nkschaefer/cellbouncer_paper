#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "Please provide num SNPs"
    exit 1
fi

nsnp=$1
bam="40k_all.bam"

{ time ./run_demux_vcf.sh $bam ${nsnp}.vcf.gz demux_vcf_donor_all_${nsnp} 2> demux_vcf_donor_all_${nsnp}.stderr; } 2> demux_vcf_donor_all_${nsnp}.time

