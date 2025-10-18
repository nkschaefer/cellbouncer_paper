## Code related to tetraploid data analysis

To do upstream:

### Create VCF file

Obtain WGS data from SRA for the following cell lines:

| line name | species |
-----------------------
| C3624 | Chimpanzee | 
| C3651 | Chimpanzee | 
| C6007B | Chimpanzee | 
| C40280 | Chimpanzee | 
| C40670 | Chimpanzee | 
| C8861 | Chimpanzee | 
| H20961 | Human | 
| H21792 | Human | 
| H23555 | Human | 
| H28126 | Human | 
| H29089 | Human | 
| CongoA4B | Bonobo |

Align and variant call (see [pipeline on GitHub](https://github.com/nkschaefer/align_pipelines)

Filter resulting VCF: 
```
bcftools view -v snps -m 2 -M 2 -c 1 -C 11 -O z -i "F_MISSING < 0.5" [variants.vcf.gz] > [variants.filtered.vcf.gz]
tabix -p vcf [variants.filtered.vcf.gz]
```

### Obtain & process RNA-seq data

Download scRNA-seq data from SRA

Run CellRanger on all libraries

Run `demux_tags` on all MULTIseq data (`multiseq` directory has well -> UID mapping)
Be sure to include library names in unique barcode identifiers
For each library:
```
demux_tags -1 [reads_R1_001.fastq.gz] -2 [reads_R2_001.fastq.gz] -N [msmap] -n [libname] -o [libname]
cat [libname].assignments >> ids_multiseq.txt 
```

Run `demux_vcf` on all 10X BAMs using VCF generated above
Again, include library name in unique IDs and concatenate all results.
Additionally, you should supply the "allowed IDs" file for the corresponding library. 
This is to prevent CellBouncer from assigning a doublet combination (tetraploid ID) that was
not included in the experiment. This is passed through the `-I` parameter and should be either
`allowed_ids_A.txt` (library A1/2) or `allowed_ids_B.txt` (library B).
For each library:
```
demux_vcf -b [possorted_genome_bam.bam] -v [variants.vcf.gz] -n [libname] -I [allowed_ids]
cat [libname].assignments >> ids_vcf.txt
```

### Extract mitochondrial & total counts from scRNA-seq matrices

This will make a list of cell barcode -> mito expr counts -> genome expr counts
But entries will only include cells with < 20% mitochondrial reads

For each library:
```
./get_mito_expr_sum.R [cellranger_outdir] [libname] >> mito_expr_sum.txt
./get_mito_expr.R [cellranger_outdir] [libname] >> mito_expr.txt
```

### Do QC and compile data

This step will load everything, filter to cells for which VCF and MULTIseq-based IDs agree, eliminate 
cells with high % mitochondrial reads, and match barcode to species, and then create necessary files
for plots.
```
./compile_data.R
```

### Infer ambient RNA & create denoised expression matrices

The above step creates `[lib]_clean.assignments` files for each library, in which low-quality cell barcodes, MULTIseq doublets, and cell barcodes with mismatching MULTIseq + VCF-based IDs are removed. These should be used for ambient RNA removal.

First, you need to copy all other relevant files for the library to match the new base file name (`[lib]_clean`):

```
cp [orig output base].counts [lib]_clean.counts
cp [orig output base].condf [lib]_clean.condf
cp [orig output base].samples [lib]_clean.samples
```

Then run `quant_contam`, passing the gene expression data.
Again, pass the correct `[allowed_ids_*.txt]` file to restrict possible identities of
tetraploid cells.

```
${CELLBOUNCER}/quant_contam -o [lib]_clean -F [CellRanger_out]/outs/filtered_feature_bc_matrix/features.tsv.gz -B [CellRanger_out]/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -M [CellRanger_out]/outs/filtered_feature_bc_matrix/matrix.mtx.gz -T [threads] -I [allowed_ids]
```

### Run DESeq2

This will assume that file names are as above: e.g. each library has a corresponding `quant_contam` run with base output name `[libname]_clean`.

Pass the file you created from the `get_mito_expr_sum.R` runs as the first argument, and then all library names
as following arguments.

```
./deseq.R [expr_sum_file] [lib1] [lib2] [lib3]
./deseq_mito_norm.R [expr_sum_file] [lib1] [lib2] [lib3]
```

### Compile results

This will create a file called `mito_expr_dat.tsv`.

```
./extract_results.R
```

### Get expression correlation of neighboring MT heavy strand genes

Where `[mito_expr]` and `[mito_expr_sum]` were created earlier

```
./corr_flanking.R [mito_expr] [mito_expr_sum]
./corr_flanking2.py
./cast_corr.R
```

### Plot stuff

```
./plot_lfc_violin.R
./plot_corr_dots.R
./plot_strand.R
./cis_trans_plot.R
```

### For comparing genotype-based demultiplexing tools

Relevant scripts for running non-CellBender genotype-based demultiplexing tools and comparing results can be found in the directory `../NSCLC`.

### For comparing ambient RNA profiling tools

Relevant scripts for running non-CellBender ambient RNA inference tools and extracting results can be found in the directory `../NSCLC/ambient/expression/`.
