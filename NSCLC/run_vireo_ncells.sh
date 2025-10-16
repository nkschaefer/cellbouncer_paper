#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "args: cells"
    exit 1
fi

cells=$1
source ~/miniforge3/bin/activate vireo
bam="40k_sc_${cells}.chosen.bam"

{ time ./run_vireo.sh $bam samp${cells}.barcodes2 500k.vcf.gz vireo_sc_${cells} 2> vireo_sc_${cells}.stderr; } 2> vireo_sc_${cells}.time

rm -r vireo_sc_${cells}_csl
rm -r vireo_sc_${cells}_vireo

