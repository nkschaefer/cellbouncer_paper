#! /usr/bin/env bash

for d in $(seq 1 7); do
    echo "${CELLBOUNCER}/demux_vcf -v ../500k.vcf.gz -b d${d}c0.05_merged.bam -o d${d}c0.05"
    echo "${CELLBOUNCER}/demux_vcf -v ../500k.vcf.gz -b d${d}c0.1_merged.bam -o d${d}c0.1"
    echo "${CELLBOUNCER}/demux_vcf -v ../500k.vcf.gz -b d${d}c0.15_merged.bam -o d${d}c0.15"
    echo "${CELLBOUNCER}/demux_vcf -v ../500k.vcf.gz -b d${d}c0.2_merged.bam -o d${d}c0.2"
    echo "${CELLBOUNCER}/demux_vcf -v ../500k.vcf.gz -b d${d}c0.25_merged.bam -o d${d}c0.25"
done
