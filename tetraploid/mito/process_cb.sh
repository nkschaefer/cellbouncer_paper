#! /usr/bin/env bash

if [ $# -lt 3 ]; then
    >&2 echo "ARGS: libname bamfile .assignments"
    exit 1
fi

lib=$1
bam=$2
assn=$3

${CELLBOUNCER}/demux_mt -b $bam -a $assn -H mito_CBfixed.haps -v mito_CBfixed.vars -i mito_CBfixed.ids -o test | bedtools intersect -wa -wb -a stdin -b mito.bed | ./sum_ase.R | sed "s/$/\t${lib}/"


