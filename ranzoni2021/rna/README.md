## Ranzoni et al 2021 RNA-seq data: for quant_contam testing

Get & prepare data
```
./prepare_data.sh
```

This will create a file of reads for Cellranger: reads_cr

Next, run CellRanger on the reads in this file (I used the 2024 human reference from the 10X website)

Then, run 
```
./after_cellranger.sh [cellranger_output_dir]
```

