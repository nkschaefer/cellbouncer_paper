#! /usr/bin/env bash

./compile_mprof.sh
./make_plot_index.R
./make_plot.R
./plot_bars.R

