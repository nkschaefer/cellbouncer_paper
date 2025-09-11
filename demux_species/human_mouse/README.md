## Code related to processing pooled human/mouse data

* Download data set from [here](https://www.10xgenomics.com/datasets/5k-hgmm-3p-nextgem)
  * Get "Gene expression - Sample barcodes" (csv file) from "Output and supplemental files" section
  * Get FASTQs from "Input files" section
    * Unpack FASTQ directory within this directory
* Download human and mouse 10X references from [here](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads)
  * The manuscript used 2020 versions of both
  * Unpack in this directory
* Download CellRanger from [here](https://www.10xgenomics.com/support/software/cell-ranger/latest)
  * Unpack anywhere
* Set up
  ```
  export CELLBOUNCER="/path/to/cellbouncer"
  ```
* Run reference building & alignment with STAR (and profile memory usage):
  ```
  ./run_starsolo_all.sh [path to human ref] [path to mouse ref] [path to cellranger]
  ```
* Run indexing & alignment with CellBouncer (and profile memory usage):
  ```
  ./run_cellbouncer_all.sh
  ```
* Plot everything
  ```
  ./plot_all.sh
  ```
  * This will make 3 plots: 
    * Memory usage for index building
    * Memory usage for species identification (and STAR alignment)
    * Proportions of cells assigned to each species using each strategy 
