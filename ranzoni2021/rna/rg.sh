#! /usr/bin/env bash
if [ $# -lt 1 ]; then
    >&2 echo "ARGS: idx"
    exit 1
fi
idx=$1

bio_misc/bin/bam_dummy_rg -b Sample${idx}_merged.bam -o Sample${idx}_merged_rg.bam -s Sample${idx} -l Sample${idx}
samtools index Sample${idx}_merged_rg.bam

