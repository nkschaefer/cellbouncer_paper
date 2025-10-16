#! /usr/bin/bash

if [ $# -lt 4 ]; then
	>&2 echo "ARGS: BAMfile bcfile VCFfile outprefix"
    exit 1 
fi

bam=$1
bc=$2
vcf=$3
outpre=$4

source /home/nkschaefer/miniforge3/etc/profile.d/conda.sh
conda activate popscle

popscle dsc-pileup --group-list $bc --sam $bam --vcf $vcf --out ${outpre}_plp
popscle freemuxlet --plp ${outpre}_plp --out ${outpre}_freemuxlet --nsample 7 

# Parse output
./freemuxlet2assn.sh ${outpre} > ${outpre}.assignments
