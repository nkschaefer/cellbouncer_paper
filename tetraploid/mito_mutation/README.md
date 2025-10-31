## Code related to the GnomAD/ClinVar plot and RNA folding

### GnomAD

Navigate to [gnomad.broadinstitute.org](https://gnomad.broadinstitute.org) and, in the search bar, type in "M-1-16559"

On the table, select "Export variants to CSV" and save as a file.

Convert CSV to TSV:

```
cat [filename].csv | sed 's/,/\t/g' > [filename].tsv
```

Make plot:
```
./gnomad_plot.R [filename].tsv
```

### RNAFold

Get human chrM fasta
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz
gunzip chrM.fa.gz
```


