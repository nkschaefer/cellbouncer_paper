#! /usr/bin/env Rscript
library(rhdf5)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("Please provide .h5 output", stderr())
    q()
}
f <- H5Fopen(args[1])

df <- data.frame(gene=f$matrix$features$name, x=f$global_latents$ambient_expression)

write.table(df, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

