#! /usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("ARGS: file.counts", stderr())
    q()
}

library(deMULTIplex2)

counts <- read.table(args[1], sep='\t', header=T)

rownames(counts) <- counts$barcode
counts <- as.matrix(counts[,-c(1)])

sink(nullfile())
res <- demultiplexTags(counts, plot.umap='none')
sink()

assn <- data.frame(barcode=rownames(res$assign_table), id=res$assign_table$final_assign)

write.table(assn, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


