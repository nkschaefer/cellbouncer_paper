## Analyses related to contamination rates & gene expression profiles of ambient RNA

Before running, first download the NSLC data raw and filtered gene expression matrices
(in MEX format, not as h5 files). For links, see two directories up (`NSCLC` directory).

Tell scripts where to find CellBouncer
```
export CELLBOUNCER="/path/to/cellbouncer"
```

Also download a human cellranger reference, and unpack it. 

Also obtain [lncAtlas](http://lncatlas.crg.eu) data in tab-separated text format, and name
it lncAtlas.txt in this directory.

Run `./run.sh`

If you run with no arguments, it will tell you the correct order of arguments.

