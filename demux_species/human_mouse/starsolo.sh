#! /usr/bin/bash

r1="5k_hgmm_3p_nextgem_fastqs/ALL_R1.fastq.gz"
r2="5k_hgmm_3p_nextgem_fastqs/ALL_R2.fastq.gz"
wl="wl.txt"
ref="mouse_human_2020"

mprof run -o mprofile_star_map.dat STAR --genomeDir $ref \
 --readFilesIn $r2 $r1 \
 --readFilesCommand zcat\
 --clipAdapterType CellRanger4\
 --outBAMsortingThreadN 1 \
 --outFileNamePrefix mouse_human_STAR \
 --outSAMattributes NH HI AS nM CR CY UR UY GX GN CB UB \
 --outSAMtype BAM SortedByCoordinate \
 --soloType CB_UMI_Simple \
 --soloCBstart 1 \
 --soloCBlen 16 \
 --soloUMIstart 17 \
 --soloUMIlen 10 \
 --soloCBwhitelist $wl \
 --outFilterScoreMin 30 \
 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
 --soloUMIfiltering MultiGeneUMI_CR \
 --soloUMIdedup 1MM_CR \
 --soloCellFilter EmptyDrops_CR \
 --soloBarcodeReadLength 0 \
 --limitSjdbInsertNsj 5000000 \
 --soloFeatures GeneFull_Ex50pAS \
 --soloMultiMappers EM 
 
