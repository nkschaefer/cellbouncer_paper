#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "args: downsample"
    exit 1
fi

ds=$1
source /home/nkschaefer/miniforge3/etc/profile.d/conda.sh
conda activate demuxalot
bam="40k_all_s${ds}.bam"

{ time ./run_demuxalot.py $bam donor_all.bcs 500k.vcf.gz 1> demuxalot_ds${ds}.assignments 2> demuxalot_ds${ds}.stderr; } 2> demuxalot_ds${ds}.time
