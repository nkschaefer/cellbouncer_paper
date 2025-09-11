#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "Please provide the path to CellRanger (unpacked)"
    exit 1
fi

cellranger=$1

./run_demux_species.sh 15 -1 $cellranger
./run_demux_species.sh 20 -1 $cellranger
./run_demux_species.sh 25 -1 $cellranger
./run_demux_species.sh 35 -1 $cellranger
./run_demux_species.sh 45 -1 $cellranger
./run_demux_species.sh 20 1000000 $cellranger
./run_demux_species.sh 20 5000000 $cellranger
./run_demux_species.sh 20 10000000 $cellranger
./run_demux_species.sh 20 20000000 $cellranger


