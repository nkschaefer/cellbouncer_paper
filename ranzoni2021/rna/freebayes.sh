#! /usr/bin/env bash
if [ $# -lt 1 ]; then
    >&2 echo "ARGS: cellranger ref"
    exit 1
fi
ref="${1}/fasta/genome.fa"
freebayes -f $ref -L bamlist.txt > vars.vcf
bgzip vars.vcf
tabix -p vcf vars.vcf.gz

bcftools view -m 2 -M 2 -v snps -i 'INFO/DP >= 100' -O z vars.vcf.gz > vars.filt.vcf.gz
tabix -p vcf vars.filt.vcf.gz


