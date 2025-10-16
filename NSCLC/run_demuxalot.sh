#! /usr/bin/env bash

if [ $# -lt 4 ]; then
    >&2 echo "USAGE: bam barcodes VCF outpre"
    exit 1
fi

source ~/miniforge3/bin/activate demuxalot
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

${SCRIPT_DIR}/run_demuxalot.py $1 $2 $3 > ${4}.assignments

