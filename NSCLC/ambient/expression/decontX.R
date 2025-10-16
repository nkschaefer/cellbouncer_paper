#! /usr/bin/env Rscript
library(Matrix)

args <- commandArgs(trailinOnly=TRUE)
if (length(args) < 1){
    write("Please provide path to filtered_feature_bc_matrix", stderr())
    q()
}
bcfile <- paste(args[1], '/barcodes.tsv.gz', sep='')
featurefile <- paste(args[1], '/features.tsv.gz', sep='')
mtxfile <- paste(args[1], 'matrix.mtx.gz', sep='')

barcodes <- read.table(bcfile)
features <- read.table(featurefile, sep='\t')
mat <- readMM(mtxfile)

# Limit to only gene expression features
rownames(mat) <- features$V2
mat <- mat[which(rownames(mat) %in% features[which(features$V3=="Gene Expression"),]$V2),]
colnames(mat) <- gsub("-1", "", barcodes$V1)

library(decontX)
sce <- decontX(mat)


contam <- data.frame(bc=colnames(mat), c=sce$contamination)
clust <- data.frame(bc=colnames(mat), clust=sce$decontX_clusters)
write.table(contam, file="decontX_contam.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(clust, file="decontX_clusters.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

saveRDS(sce$decontXcounts, file="decontX_counts_adj.RDS")

# https://github.com/campbio/decontX/discussions/12
eta <- metadata(sce)$estimates$all_cells$eta
rownames(eta) <- rownames(sce$decontXcounts)
contam_prof <- data.frame(gene=rownames(eta), x=rowSums(eta))
contam_prof$x <- contam_prof$x / sum(contam_prof$x)

write.table(contam_prof, file="decontX_genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

