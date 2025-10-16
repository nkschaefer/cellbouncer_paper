#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "ERROR: please provide .counts"
    exit 1
fi

if [ $( echo -e "$1" | grep "counts" | wc -l ) -eq 0 ]; then
    >&2 echo "ERROR: $1 is not a .counts file"
    exit 1
fi

counts=$( echo -e "${1%.counts}" | sed 's/\./_/' )
counts=${counts}.counts
base=${counts%.counts}
if [ "$1" != "$counts" ]; then
    cp $1 $counts
fi

source /Users/nathan/opt/miniconda3/etc/profile.d/conda.sh
conda activate multiseq

# Creates mtx format stuff
./run_GMMdemux.sh $counts | ./hetdoub.py | sed "s/^/${base}\tGMMdemux\t/"
./run_demuxmix.R $counts | ./hetdoub.py | sed "s/^/${base}\tdemuxmix\t/"
./run_BFF.R $counts | ./hetdoub.py | sed "s/^/${base}\tBFF\t/"
./run_hashedDrops.R $counts | ./hetdoub.py | sed "s/^/${base}\thashedDrops\t/"

conda deactivate
conda activate signac
./run_HTOdemux.sh $counts | ./hetdoub.py | sed "s/^/${base}\tHTOdemux\t/"

conda deactivate
./run_demultiplex2.R $counts | ./hetdoub.py | sed "s/^/${base}\tdeMULTIplex2\t/"

conda activate solo
./run_hashsolo.py ${counts%.counts}_mtx | ./hetdoub.py | sed "s/^/${base}\thashSolo\t/"

./refmt_assign.sh ${counts%.counts}.assignments | ./hetdoub.py | sed "s/^/${base}\tCellBouncer\t/"


