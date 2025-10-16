#! /usr/bin/env bash

if [ "$1" == "" ]; then
    >&2 echo "ERROR: please provide cells"
    exit 1
fi

cells=$1

source /home/nkschaefer/miniforge3/etc/profile.d/conda.sh
conda activate demuxalot

bam="40k_sc_${cells}.chosen.bam"

{ time ./run_demuxalot.py $bam samp${cells}.barcodes2 500k.vcf.gz 1> demuxalot_sc_${cells}.assignments 2> demuxalot_sc_${cells}.stderr; } 2> demuxalot_sc_${cells}.time

