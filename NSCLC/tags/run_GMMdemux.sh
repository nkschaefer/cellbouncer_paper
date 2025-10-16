#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "USAGE: run_GMMdemux.sh .counts"
    exit 1
fi

countsfn=$1

if [ -d ${countsfn%.counts}_mtx ]; then
    rm -r ${countsfn%.counts}_mtx
fi

./counts2mtx.R $countsfn > /dev/null
gzip ${countsfn%.counts}_mtx/features.tsv
gzip ${countsfn%.counts}_mtx/barcodes.tsv
gzip ${countsfn%.counts}_mtx/matrix.mtx

cols=$( cat $countsfn | head -1 | cut -f2- | sed 's/\t/,/g' )

GMM-demux ${countsfn%.counts}_mtx $cols -x $cols -f ${countsfn%.counts}_gmmdemux > /dev/null

./parse_GMMdemux.R ${countsfn%.counts}_gmmdemux


