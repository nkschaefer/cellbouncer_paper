#! /usr/bin/env bash
if [ $# -lt 1 ]; then
    >&2 echo "ARGS: index"
    exit 1
fi

idx=$1

cd Sample${idx}_bc
outbase="Ranzoni_S${idx}_L001"
r1="${outbase}_R1_001.fastq"
r2="${outbase}_R2_001.fastq"
if [ -e "${r1}" ]; then
    rm $r1
fi
if [ -e "${r2}" ]; then
    rm $r2
fi
if [ -e "${r1}.gz" ]; then
    rm "${r1}.gz"
fi
if [ -e "${r2}.gz" ]; then
    rm "${r2}.gz"
fi
ls *_R1*.fastq.gz | sort -V | while read fn; do
    fn2=$( echo -e "$fn" | sed 's/_R1/_R2/' )
    zcat $fn >> $r1
    zcat $fn2 >> $r2
done

gzip $r1
gzip $r2

