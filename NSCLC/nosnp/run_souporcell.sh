#! /usr/bin/env bash

if [ $# -lt 3 ]; then
    >&2 echo "USAGE: run_souporcell.sh bam barcodes out"
    exit 1
fi

bam=$1
barcodes=$2
outpre=$3

source /home/nkschaefer/miniforge3/etc/profile.d/conda.sh

conda activate souporcell
export PATH="${PATH}:."  
souporcell/souporcell_pipeline.py -i $bam -b $barcodes -f hg38.fa -t 8 -o $outpre -k 7
souporcell2assn.sh ${outpre} > ${outpre}.assignments

