#! /usr/bin/env bash

# Creates a directory with symlinks in the format expected by the demux_species pipeline
if [ ! -d reads ]; then
    mkdir reads
fi
wd=$(pwd)

find . -mindepth 1 -maxdepth 1 -type d -name "SAMN*" | while read sample; do
    sample="${sample##*/}"
    si=1
    ls ${sample}/*_1.fastq.gz | while read fn; do
        fn2="${fn##*/}"
        prefix="${fn2%%_*.fastq.gz}"
        suffix="${fn2##*_}"
        suffix="${suffix%%.fastq.gz}"
        fn_new="reads/${sample}_S${si}_L00${si}_R1_001.fastq.gz"
        fn_newr2="reads/${sample}_S${si}_L00${si}_R2_001.fastq.gz"
        ln -s ${wd}/$fn $fn_new
        fnr2=$( echo -e "${fn}" | sed 's/_1/_2/' )
        ln -s ${wd}/$fnr2 $fn_newr2
        si=$(( $si + 1 ))
    done
done
