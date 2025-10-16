#! /usr/bin/bash
if [ $# -lt 1 ]; then
    >&2 echo "ARGS: idx"
    exit 1
fi
idx=$1
samp="Sample${idx}"

sampfile="${samp}_files.txt"
ls mapped/${samp}/*.bam > $sampfile

bio_misc/merge_many_bams.sh $sampfile "${samp}_merged.bam"
samtools index ${samp}_merged.bam

