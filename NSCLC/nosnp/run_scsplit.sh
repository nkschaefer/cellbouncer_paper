#! /usr/bin/env bash

if [ $# -lt 4 ]; then
    >&2 echo "ARGS: BAM VCF barcode out"
    exit 1
fi

bam=$1
vcf=$2
bc=$3
outpre=$4

source /home/nkschaefer/miniforge3/etc/profile.d/conda.sh
conda activate scSplit

if [ ! -d $outpre ]; then
    mkdir $outpre
fi

time scSplit/scSplit count -v $vcf -i $bam -b $bc -r ref.tsv -a alt.tsv -o $outpre

# Run scSplit
scSplit/scSplit run -r ${outpre}/ref.tsv -a ${outpre}/alt.tsv -n 7 -o ${outpre}

# Reformat output
./scSplit2assn.R ${outpre} > ${outpre}.assignments

