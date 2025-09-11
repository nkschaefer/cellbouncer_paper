### Code for processing data from [Edenhofer et al (2024)](https://pubmed.ncbi.nlm.nih.gov/38947524/)

These data consist of scRNA-seq from human and macaque cells, identifiable by the CRISPR sgRNAs targeted to the cells.

Some guides targeted both species and these could not be used to indicate species of origin.

This analysis expects data to be located in a subdirectory called `geo_submission_hellmann_18082023`, containing reads (in FASTQ format) as well as the file `protospacer_calls_per_cell.csv`. Download the relevant file and unpack in this directory.

Output: 

* a PNG and PDF plot showing concordance of sgRNA-based IDs (Y) and demux_species k-mer-based IDs (X)

* a file (precision_recall.txt) storing the precision and recall of demux_species assignments treating the sgRNA-based assignments as truth

#### To run
```
export CELLBOUNCER="/path/to/cellbouncer"
./make_index.sh
./run_demux_species.sh
```
