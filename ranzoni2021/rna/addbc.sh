#! /usr/bin/env bash
if [ $# -lt 1 ]; then
    >&2 echo "ARGS: index"
    exit 1
fi
idx=$1

uid=$( head -n $idx rnaseq_ids_bc_cellranger.txt | tail -1 | cut -f1 )
bc=$( head -n $idx rnaseq_ids_bc_cellranger.txt | tail -1 | cut -f2 )
samp=$( grep -w $uid rnaseq_files.txt | tail -1 | cut -f2 )
run=$( grep -w $uid samp2run.txt | head -1 | cut -f2 )

if [ ! -d "${samp}_bc" ]; then
   mkdir "${samp}_bc"
 fi

./addbc.py -1 "${samp}/${run}_1.fastq.gz" -2 "${samp}/${run}_2.fastq.gz" -b $bc -o "${samp}_bc/${run}"  
