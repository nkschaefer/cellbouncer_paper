#! /usr/bin/env Rscript
library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1){
    write("Please provide cellranger output dir", stderr())
    q()
}
dirn <- args[1]

barcodes <- read.table(paste(dirn, '/outs/raw_feature_bc_matrix/barcodes.tsv.gz', sep=''))
features <- read.table(paste(dirn, '/outs/raw_feature_bc_matrix/features.tsv.gz', sep=''), sep='\t')
mat <- readMM(paste(dirn, '/outs/raw_feature_bc_matrix/matrix.mtx.gz', sep=''))

rownames(mat) <- features$V2
mat <- mat[which(rownames(mat) %in% features[which(features$V3=="Gene Expression"),]$V2),]
colnames(mat) <- gsub("-1", "", barcodes$V1)

#barcodesRaw <- read.table('outs/raw_feature_bc_matrix/barcodes.tsv.gz')
#featuresRaw <- read.table('outs/raw_feature_bc_matrix/features.tsv.gz', sep='\t')
#matRaw <- readMM('outs/raw_feature_bc_matrix/matrix.mtx.gz')

#rownames(matRaw) <- featuresRaw$V2
#matRaw <- matRaw[which(rownames(matRaw) %in% featuresRaw[which(featuresRaw$V3=="Gene Expression"),]$V2),]
#colnames(matRaw) <- gsub("-1", "", barcodesRaw$V1)

library(decontX)
sce <- decontX(mat)

diff <- mat - sce$decontXcounts
diffGenes <- data.frame(gene=rownames(mat), diff=rowSums(diff))
diffGenes$orig <- rowSums(mat)


contam <- data.frame(bc=colnames(mat), c=sce$contamination)

write.table(contam, file="decontX_contam.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#write.table(diffGenes, file="decontX_genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

saveRDS(sce$decontXcounts, file="decontX_counts_adj.RDS")

clusts <- data.frame(bc=colnames(mat), clust=sce$z)
write.table(clusts, file="decontX_clusts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# https://github.com/campbio/decontX/discussions/12
eta <- sce$estimates$all_cells$eta
rownames(eta) <- rownames(sce$decontXcounts)
contam_prof <- data.frame(gene=rownames(eta), x=rowSums(eta))
contam_prof$x <- contam_prof$x / sum(contam_prof$x)

write.table(contam_prof, file="decontX_contam_prof.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

