#! /usr/bin/env bash

if [ $# -lt 3 ]; then
    >&2 echo "ARGS: bam vcf out"
    exit 1
fi

bam=$1
vcf=$2
outpre=$3
${CELLBOUNCER}/demux_vcf -b $bam -o $outpre -v $vcf

