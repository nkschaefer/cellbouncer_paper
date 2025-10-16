#! /usr/bin/env bash

zcat Sample1_bc/*_1.fastq.gz | gzip -c - > Sample1_R1.fq.gz
zcat Sample1_bc/*_2.fastq.gz | gzip -c - > Sample1_R2.fq.gz

r1="Sample1_R1.fq.gz"
r2="Sample1_R2.fq.gz"
ref="hg38.idx"
lib="Sample1"

minimap2 -y -a -x sr -R "@RG\tID:${lib}\tSM:${lib}\tPL:Illumina" $ref $r1 $r2 | samtools sort -o ${lib}_sorted.bam
samtools index ${lib}_sorted.bam

