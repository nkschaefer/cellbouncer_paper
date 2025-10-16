#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "args: downsample"
    exit 1
fi

ds=$1

source ~/miniforge3/bin/activate vireo
bam="40k_all_s${ds}.bam"

{ time ./run_vireo.sh $bam donor_all.bcs 500k.vcf.gz vireo_ds${ds} 2> vireo_ds${ds}.stderr; } 2> vireo_ds${ds}.time

rm -r vireo_ds${ds}_csl
rm -r vireo_ds${ds}_vireo
