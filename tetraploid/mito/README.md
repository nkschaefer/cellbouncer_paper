## Code related to mitochondrial identification

To cluster (and subcluster) mitochondrial haplotypes:
Run `${CELLBOUNCER}/demux_mt`, followed by `${CELLBOUNCER}/mt_subcluster.py` the number of times specified
in the manuscript.

Plots comparing "true" mitochondrial haplotypes (from WGS variants) and the results of clustering were made using the script `${CELLBOUNCER}/plot/compare_assignments.R`.

### VCF-based mitochondrial haplotype assignment

Three files were used for assigning mitochondrial haplotypes, using variants learned
from WGS + variant calling:

```
mito.vars
mito.haps
mito.ids
```

These can be fed to demux_mt: 
```
demux_mt -b [bamfile] -v mito.vars -H mito.haps -i mito.ids
```

To count reads matching each mitochondrial haplotype in each cell, run this along with
the cell-ID mapping learned from demux_vcf:
```
demux_mt -b [bamfile] -a [demux_vcf_assignments] -v mito.vars -H mito.haps -i mito.ids > [mito_counts]
```

This will give you a file mapping cell barcode -> reads matching mitochondrial haplotype 1 -> reads matching mitochondrial haplotype 2

The above will give you counts of reads matching one or the other specific haplotype.

For the sake of the paper, we used species-fixed differences in mitochondrial haplotypes to determine
presence of species' mitochondrial haplotypes. The relevant files for human vs chimp haplotypes are:

```
mito_HCfixed.vars
mito_HCfixed.haps
mito_HCfixed.ids
```

And for chimp vs bonobo:
```
mito_CBfixed.vars
mito_CBfixed.haps
mito_CBfixed.ids
```

After running demux_vcf and demux_tags on our data, we chose cells with agreement between the VCF-based
and tag-based ID, and then split cells by species of origin (human/human, chimp/chimp, human/chimp, chimp/bonobo)
and then counted mitochondrial reads in each separately:

for human/human: use `mito.*` files
for chimp/chimp: use `mito.*` files
for human/chimp: use `mito_HCfixed.*` files
for chimp/bonobo: use `mito_CBfixed.*` files


The above can be done by these scripts, which can be used to compile data:
**NOTE**: the `[assignments]` files listed below come from `demux_vcf`, but they should also 
be filtered to contain only barcodes for which the `demux_vcf` assignments (from genotypes)
match the `demux_tags` assignments (from MULTIseq data). See parent directory.

If you have already run code in the parent directory (including the `compile_assignments.R` script), there should be files called `[libname]_hh.assignments`, `[libname]_hc.assignments`, `[libname]_cc.assignments`, and `[libname]_cb.assignments`, where `hh` is human/human, `hc` is human/chimp, `cc` is chimp/chimp, and `cb` is chimp/bonobo.

```
# For human/human fusions:
./process_auto.sh [libname] [bam] [assignments] >> hh_data.txt
# For chimp/chimp fusions:
./process_auto.sh [libname] [bam] [assignments] >> cc_data.txt
# For human/chimp fusions:
./process_hc.sh [libname] [bam] [assignments] >> hc_data.txt
# For chimp/bonobo fusions:
./process_cb.sh [libname] [bam] [assignments] >> cb_data.txt
```

Next, compile data:
```
./compile_data.R
```

This will create `inter_species_dat.tsv` and `intra_species_dat.tsv`


Then, use the resulting count data to assign species of origin using mixture models:
**NOTE**: this requires an old version of pomegranate, before API changes (I used 0.15.0)

```
./assign_mito.py
```

This will create `mito_ase_results_intra_species.tsv`, `mito_ase_results_HC.tsv`, and `mito_ase_results_CB.tsv`.

### Create figures

Plot showing proportions of mito1/(mito1+mito2) in each cell
```
./plot_mito_frac.R
```

Reformat mito haplotypes to make them easier to load
```
cat mito.haps | sed 's/./&\t/g > mito.haps2
```

Make Fig S13 plot
```
./make_figs13.R 
```
 
