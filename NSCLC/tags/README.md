## Code related to testing tag (in this case, 10X CellPlex) processing

Before running, obtain reads from 10X NSCLC experiment (see above directory)
and unpack here.

This should create these directories, within this directory:
```
40k_NSCLC_DTC_3p_HT_nextgem_cmo_1_fastqs
40k_NSCLC_DTC_3p_HT_nextgem_cmo_2_fastqs
40k_NSCLC_DTC_3p_HT_nextgem_cmo_3_fastqs
```

Tell scripts where to find CellBouncer
```
export CELLBOUNCER="/path/to/cellbouncer"
```

Then:
```
./run_everything.sh [/path/to/cellranger/barcode_whitelist]
```

