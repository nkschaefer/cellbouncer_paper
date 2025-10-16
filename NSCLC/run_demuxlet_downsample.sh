#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "args: downsample"
    exit 1
fi

ds=$1
source ~/miniforge3/bin/activate popscle
bam="40k_all_s${ds}.bam"

{ time ./run_demuxlet.sh $bam 500k.vcf.gz demuxlet_ds${ds} 2> demuxlet_ds${ds}.stderr; } 2> demuxlet_ds${ds}.time

