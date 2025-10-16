#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "Please provide num SNPs"
    exit 1
fi
nsnp=$1
source ~/miniforge3/bin/activate vireo
bam="40k_all.bam"

{ time ./run_vireo.sh $bam donor_all.bcs ${nsnp}.vcf.gz vireo_donor_all_${nsnp} 2> vireo_donor_all_${nsnp}.stderr; } 2> vireo_donor_all_${nsnp}.time

rm -r vireo_donor_all_${nsnp}_csl
rm -r vireo_donor_all_${nsnp}_vireo
