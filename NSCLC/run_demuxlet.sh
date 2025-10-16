#! /usr/bin/env bash

if [ $# -lt 3 ]; then
	>&2 echo "ARGS: bam vcf out"
	exit 1
fi

bam=$1
vcf=$2
outpre=$3

source ~/miniforge3/bin/activate popscle

popscle demuxlet --vcf $vcf --field GT --sam $bam --out $outpre

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

${SCRIPT_DIR}/demuxlet2assn.sh ${outpre}.best > ${outpre}.assignments

 
