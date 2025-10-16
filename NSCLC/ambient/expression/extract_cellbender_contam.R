#! /usr/bin/env Rscript

library(rhdf5)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("ERROR: please provide output (filtered) h5", stderr())
    q()
}

h5 <- args[1]
f <- H5Fopen(h5)
percs <- data.frame(bc=f$matrix$barcodes, c=f$droplet_latents$background_fraction)

write.table(percs, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

