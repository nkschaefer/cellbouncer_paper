#! /usr/bin/bash

if [ $# -lt 2 ]; then
    >&2 echo "args: cellranger_ref index"
    exit 1
fi
ref=$1
idx=$2
sample="Sample${idx}"

grep -w $sample ../rnaseq_files.txt | cut -f1 | uniq | while read id; do
    
    bc=$( grep -w $id rnaseq_ids_bc.txt | head -1 | cut -f2 )
    f1=$( grep -w $id rnaseq_files.txt | head -1 | cut -f3 )
    f2=$( grep -w $id rnaseq_files.txt | tail -1 | cut -f3 )

    f1="${f1##*/}"
    f2="${f2##*/}"
    
    samp=$sample

    if [ ! -d "mapped/${samp}" ]; then
        mkdir "mapped/${samp}"
    fi
   
    ref="${ref}/star/" 
    star="STAR"
    
    if [ ! -e "mapped/${samp}/${id}.bc.bam" ]; then
        $star --runThreadN 1 \
            --outFileNamePrefix "mapped/${samp}/${id}" \
            --outSAMtype BAM SortedByCoordinate \
            --genomeDir $ref \
            --readFilesIn "${samp}/${f1}" "${samp}/${f2}" \
            --readFilesCommand zcat
         
        rm "mapped/${samp}/${id}Log.final.out"
        rm "mapped/${samp}/${id}Log.progress.out"
        rm "mapped/${samp}/${id}Log.out"
        rm "mapped/${samp}/${id}J.out.tab"
        
        # Insert cell barcode
        bio_misc/bin/bam_add_tag -b "mapped/${samp}/${id}Aligned.sortedByCoord.out.bam" \
            -o "mapped/${samp}/${id}.bc.bam" \
            -t CB \
            -T $bc

        rm "mapped/${samp}/${id}Aligned.sortedByCoord.out.bam"
    fi 
    samtools index "mapped/${samp}/${id}.bc.bam"

done
