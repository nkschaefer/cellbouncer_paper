#! /usr/bin/bash

if [ $# -lt 4 ]; then
	>&2 echo "ARGS: BAMfile bcfile VCFfile outprefix"
fi

bam=$1
bc=$2
vcf=$3
outpre=$4

source ~/miniforge3/bin/activate vireo
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cellsnp-lite --gzip -s $bam -T $vcf -O ${outpre}_csl --cellTAG CB -b $bc --genotype
vireo -t GT -d $vcf -c ${outpre}_csl -o ${outpre}_vireo
${SCRIPT_DIR}/vireo2assn.sh ${outpre}_vireo/donor_ids.tsv > ${outpre}.assignments

