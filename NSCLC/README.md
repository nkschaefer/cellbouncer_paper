## Scripts for processing 10X 40k NSCLC data

Get data from here: https://www.10xgenomics.com/datasets/40-k-mixture-of-nsclc-dt-cs-from-7-donors-3-ht-v-3-1-3-1-high-6-1-0

Tell scripts where to find CellBouncer
```
export CELLBOUNCER="/path/to/cellbouncer"
```

Run `./get_data.sh` to merge BAMs, call variants, and downsample data as needed to run tests

### Genotype-based demultiplexing trials

There are three ways the data were downsampled: 
* Number of SNPs
  * Processing scripts: `run_*_nsnp.sh` 
* Reads per cell
  * Processing scripts: `run_*_downsample.sh`
* Number of cells
  * Processing scripts: `run_*_ncells.sh`

To run everything in parallel: 
```
./make_vcf_jobs.sh | parallel -j [njobs]
```
Where `[njobs]` is the number of parallel processes to run

Then create consensus IDs:
```
./get_consensus.R > ids_true_consensus.txt
```

Compile data about accuracy:
```
./sens_spec_all.sh
```
Compile data about execution time:
```
./time_all.sh
```

Make plots
```
./make_plots.R
```
