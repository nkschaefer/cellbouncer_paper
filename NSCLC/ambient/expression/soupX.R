#! /usr/bin/env Rscript
library(SoupX)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2){
    write("ARGS: filtered_feature_bc_matrix raw_feature_bc_matrix", stderr())
    q()
}

toc <- Seurat::Read10X(args[1])
tod <- Seurat::Read10X(args[2])
toc2 <- toc[['Gene Expression']]
tod2 <- tod[['Gene Expression']]

clust <- read.table('decontX_clusters.txt', sep='\t', header=T)
colnames(clust)[2] <- "ID"
rownames(clust) <- clust$bc
clust <- clust[,c(2),drop=FALSE]
assn <- clust

toc2 <- toc2[,which(colnames(toc2) %in% rownames(assn))]
assn <- assn[match(colnames(toc2), rownames(assn)),,drop=FALSE]
a2 <- setNames(assn$ID, rownames(assn))

sc <- SoupChannel(tod2, toc2)
sc <- setClusters(sc, a2)

sc <- autoEstCont(sc)

# Cell contamination rates
contam <- sc$metaData
contam$bc <- gsub("-1", "", rownames(contam))
contam <- contam[,c(4,3)]

library(Matrix)
counts_out <- adjustCounts(sc, roundToInt=TRUE)
counts_adj <- data.frame(bc=gsub("-1", "", colnames(sc$toc)), before=Matrix::colSums(sc$toc))
counts_adj$after <- Matrix::colSums(counts_out)
counts_adj$c <- (1-counts_adj$after/counts_adj$before)
contam <- merge(contam, counts_adj)

write.table(contam, file="soupX_contam.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

prof <- sc$soupProfile
prof$gene <- rownames(prof)
prof <- prof[,c(3,1,2)]
write.table(prof, file="soupX_genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

saveRDS(counts_out, file="soupX_counts_adj_dxclust.RDS")

