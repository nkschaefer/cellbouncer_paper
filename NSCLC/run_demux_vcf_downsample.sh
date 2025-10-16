#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "args: downsample"
    exit 1
fi

ds=$1

bam="40k_all_s${ds}.bam"

{ time ./run_demux_vcf.sh $bam 500k.vcf.gz demux_vcf_ds${ds} 2> demux_vcf_ds${ds}.stderr; } 2> demux_vcf_ds${ds}.time

