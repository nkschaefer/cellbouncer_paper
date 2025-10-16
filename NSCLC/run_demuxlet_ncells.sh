#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "args: cells"
    exit 1
fi

cells=$1
source ~/miniforge3/bin/activate popscle
bam="40k_sc_${cells}.chosen.bam"

{ time ./run_demuxlet.sh $bam 500k.vcf.gz demuxlet_sc_${cells} 2> demuxlet_sc_${cells}.stderr; } 2> demuxlet_sc_${cells}.time

