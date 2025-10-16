#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "USAGE: run_HTOdemux.sh .counts"
    exit 1
fi

countsfn=$1

if [ ! -d ${countsfn%.counts}_mtx ]; then
    ./counts2mtx.R $countsfn > /dev/null
fi

./run_HTOdemux.R ${countsfn%.counts}_mtx

