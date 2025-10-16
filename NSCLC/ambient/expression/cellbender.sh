#! /usr/bin/bash

if [ $# -lt 1 ]; then
    >&2 echo "Please provide path to (unpacked) raw mtx format data for NSCLC experiment"
    exit 1
fi

cellbender remove-background --input $1 \
 --output NSCLC_cellbender.h5 \
 --expected-cells 33000 \
 --model full \
 --epochs 150 \
 --cuda

# Extract relevant information
./extract_cellbender_contam.R NSCLC_cellbender.h5 > cellbender_contam.txt
./extract_cellbender_profile.R NSCLC_cellbender.h5 > cellbender_genes.txt

