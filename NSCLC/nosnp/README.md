## Demultiplexing NSLCLC data without prior knowledge of SNPs

Download freemuxlet, scSplit, souporcell, and vireo
Place/unpack scSplit and souporcell in this directory

Edit scripts in this directory to change conda paths to wherever you
have installed conda (and whatever conda env names you install these
programs' environments under)

### CellBouncer:
```
${CELLBOUNCER}/demux_mt -b ../40k_all.bam -o demux_mt -c 
```
To generate the heatmap plot:
```
${CELLBOUNCER}/plot/demux_mt demux_mt
```
### Run everything else
```
./run_everything.sh
```

### Generate plots
Need to manually compile output from the `${CELLBOUNCER}/plot/compare_assignments.R` runs (see end of `run_everything.sh` to capture the F1 statistics (if file 1 is truth)

Also need to compute number of variants in each VCF using `count_variable.py`, e.g.:
```
zcat [file.vcf.gz] | ./count_variable.py
```
Place these in a file called compare_dat.txt (check current file for format)
Run plot script:
```
./plot_compare_dat.R
```
 
