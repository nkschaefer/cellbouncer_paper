#! /usr/bin/env bash

zcat Sample2_bc/*_1.fastq.gz | gzip -c - > Sample2_R1.fq.gz
zcat Sample2_bc/*_2.fastq.gz | gzip -c - > Sample2_R2.fq.gz

r1="Sample2_R1.fq.gz"
r2="Sample2_R2.fq.gz"
ref="hg38.idx"
lib="Sample2"

minimap2 -y -a -x sr -R "@RG\tID:${lib}\tSM:${lib}\tPL:Illumina" $ref $r1 $r2 | samtools sort -o ${lib}_sorted.bam
samtools index ${lib}_sorted.bam

