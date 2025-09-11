#! /usr/bin/env bash

human_tx="gencode.v49.transcripts.fa"
mouse_tx="gencode.vM38.transcripts.fa"

k=$1

mprof run -o mprofile_index_k${k}.dat ${CELLBOUNCER}/utils/demux_species_ref.py -k $k -n Human Mouse -f $human_tx $mouse_tx -O -o hm${k} -N 20000000

