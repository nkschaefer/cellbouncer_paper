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
base2=$( echo -e "$base" | sed 's/\./_/' )
cp $1 ${base2}.counts
counts2=${base2}.counts

source /Users/nathan/opt/miniconda3/etc/profile.d/conda.sh
conda activate multiseq

# Creates mtx format stuff
./run_GMMdemux.sh $counts2 | ./assess.R | sed "s/^/${base}\tGMMdemux\t/"
./run_demuxmix.R $counts2 | ./assess.R | sed "s/^/${base}\tdemuxmix\t/"
./run_BFF.R $counts2 | ./assess.R | sed "s/^/${base}\tBFF\t/"
./run_hashedDrops.R $counts2 | ./assess.R | sed "s/^/${base}\thashedDrops\t/"

conda deactivate
conda activate signac
./run_HTOdemux.sh $counts2 | ./assess.R | sed "s/^/${base}\tHTOdemux\t/"

conda deactivate
./run_demultiplex2.R $counts2 | ./assess.R | sed "s/^/${base}\tdeMULTIplex2\t/"

conda activate solo
./run_hashsolo.py ${counts2%.counts}_mtx | ./assess.R | sed "s/^/${base}\thashSolo\t/"

${CELLBOUNCER}/demux_tags -o ${base2}
./refmt_assign.sh ${base2}.assignments | ./assess.R | sed "s/^/${base}\tCellBouncer\t/"


