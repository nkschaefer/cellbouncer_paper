#! /usr/bin/bash

if [ $# -lt 4 ]; then
	>&2 echo "ARGS: BAMfile bcfile VCFfile outprefix"
    exit 1 
fi

bam=$1
bc=$2
vcf=$3
outpre=$4

cellsnp-lite --gzip -s $bam -T $vcf -O ${outpre}_csl --cellTAG CB -b $bc --genotype
vireo -N 7 -c ${outpre}_csl -o ${outpre}_vireo
./vireo2assn.sh ${outpre}_vireo/donor_ids.tsv > ${outpre}.assignments

