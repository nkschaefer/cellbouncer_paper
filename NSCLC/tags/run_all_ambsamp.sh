#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "ERROR: please provide .counts"
    exit 1
fi

if [ $( echo -e "$1" | grep "counts" | wc -l ) -eq 0 ]; then
    >&2 echo "ERROR: $1 is not a .counts file"
    exit 1
fi

base="${1%.counts}"
counts=$1

source /Users/nathan/opt/miniconda3/etc/profile.d/conda.sh
conda activate multiseq

# Creates mtx format stuff
./run_GMMdemux.sh $counts | ./assess.R | sed "s/^/${base}\tGMMdemux\t/"
./run_demuxmix.R $counts | ./assess.R | sed "s/^/${base}\tdemuxmix\t/"
./run_BFF.R $counts | ./assess.R | sed "s/^/${base}\tBFF\t/"
./run_hashedDrops.R $counts | ./assess.R | sed "s/^/${base}\thashedDrops\t/"

conda deactivate
conda activate signac
./run_HTOdemux.sh $counts | ./assess.R | sed "s/^/${base}\tHTOdemux\t/"

conda deactivate
./run_demultiplex2.R $counts | ./assess.R | sed "s/^/${base}\tdeMULTIplex2\t/"

conda activate solo
./run_hashsolo.py ${counts%.counts}_mtx | ./assess.R | sed "s/^/${base}\thashSolo\t/"

${CELLBOUNCER}/demux_tags -o ${base}
./refmt_assign.sh ${base}.assignments | ./assess.R | sed "s/^/${base}\tCellBouncer\t/"

