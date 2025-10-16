#! /usr/bin/env bash
if [ $# -lt 2 ]; then
    >&2 echo "ARGS: sampleID index"
    exit 1
fi
sample=$1
idx=$2

r1=$( ls ${sample}/*_1.fastq.gz | sort -V | head -n $idx | tail -1 )
r2=$( echo -e "$r1" | sed 's/_1.fastq.gz/_2.fastq.gz/' )
bc=$( head -n $idx ${sample}_barcodes.txt | tail -1 )
outd="${sample}_bc"

./addbc.py $r1 $r2 $outd $bc
