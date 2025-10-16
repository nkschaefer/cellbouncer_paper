#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1){
    write("ARGS: file.counts", stderr())
    q()
}
library(Matrix)

counts <- read.table(args[1], sep="\t", header=T)

basen <- strsplit(args[1], '.', fixed=T)[[1]][1]

dir_out <- paste(basen, '_mtx', sep='')
dir.create(dir_out)

out_bc <- paste(dir_out, '/barcodes.tsv', sep='')
write.table(counts[,c(1),drop=FALSE], file=out_bc, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
out_features <- paste(dir_out, '/features.tsv', sep='')
write.table(data.frame(feature=colnames(counts[,-c(1)])), file=out_features, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

out_mtx <- paste(dir_out, '/matrix.mtx', sep='')
rownames(counts) <- counts$barcode
m <- t(as.matrix(counts[,-c(1)]))
m <- as(m, 'CsparseMatrix')
writeMM(m, out_mtx)



