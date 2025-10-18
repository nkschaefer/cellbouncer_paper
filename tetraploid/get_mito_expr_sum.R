#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2){
    write("ERROR: provide cellranger out dir and lib", stderr())
    q()
}

lib <- args[2]

library(Matrix)

bc <- read.table(paste(args[1], '/outs/filtered_feature_bc_matrix/barcodes.tsv.gz', sep=''))
gene <- read.table(paste(args[1], '/outs/filtered_feature_bc_matrix/features.tsv.gz', sep=''), sep='\t')
mat <- readMM(paste(args[1], '/outs/filtered_feature_bc_matrix/matrix.mtx.gz', sep=''))

colnames(mat) <- bc$V1
rownames(mat) <- gene$V2

mito <- mat[grep("^MT-", rownames(mat)),]

cellsums <- data.frame(bc=colnames(mat), tot=colSums(mat))
mitosums <- data.frame(bc=colnames(mito), mito=colSums(mito))

merged <- merge(mitosums, cellsums)
merged$bc <- gsub("-1", "", merged$bc)
merged$bc <- paste(merged$bc, lib, sep='-')

merged <- merged[which(merged$mito/merged$tot < 0.2),]

write.table(merged, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

