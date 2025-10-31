#! /usr/bin/env bash

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz
gunzip chrM.fa.gz

echo "Human allele"

echo -e "chrM\t10458\t10492" | bedtools getfasta -fi chrM.fa -bed stdin | rnafold -p
./extr_dot_pos.R chrM_10458-10492_dp.ps 20

echo "Ancestral allele"
begin=$( echo -e "chrM\t10458\t10492" | bedtools getfasta -fi chrM.fa -bed stdin | tail -1 | cut -c1-19 )
end=$( echo -e "chrM\t10458\t10492" | bedtools getfasta -fi chrM.fa -bed stdin | tail -1 | cut -c21-34 )
echo -e ">chrM_10458-10492_anc\n${begin}T${end}" | rnafold -p
./extr_dot_pos.R chrM_10458-10492_anc_dp.ps 20


