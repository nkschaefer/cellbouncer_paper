#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "args: cells"
    exit 1
fi

cells=$1

bam="40k_sc_${cells}.chosen.bam"

{ time ./run_demux_vcf.sh $bam 500k.vcf.gz demux_vcf_sc_${cells} 2> demux_vcf_sc_${cells}.stderr; } 2> demux_vcf_sc_${cells}.time

