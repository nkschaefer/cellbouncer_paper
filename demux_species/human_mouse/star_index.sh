#! /bin/bash

genome_fa=mouse_human_2020.fa
genome_gtf=mouse_human_2020.gtf
genome_out=mouse_human_2020

# limitRAM = memGB * 1024 ** 3
# 214748364800 is 200 GB

star="STAR"
mprof run -o mprofile_star_index.dat $star --runMode genomeGenerate \
    --sjdbGTFfile $genome_gtf \
    --genomeFastaFiles $genome_fa \
    --genomeDir $genome_out \
    --limitGenomeGenerateRAM 214748364800

