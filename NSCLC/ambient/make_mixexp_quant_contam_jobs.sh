#! /usr/bin/env bash

for d in $(seq 1 7); do
    echo "${CELLBOUNCER}/quant_contam -D 0.1 -o d${d}c0.05"
    echo "${CELLBOUNCER}/quant_contam -D 0.1 -o d${d}c0.1"
    echo "${CELLBOUNCER}/quant_contam -D 0.1 -o d${d}c0.15"
    echo "${CELLBOUNCER}/quant_contam -D 0.1 -o d${d}c0.2"
    echo "${CELLBOUNCER}/quant_contam -D 0.1 -o d${d}c0.25"
done
