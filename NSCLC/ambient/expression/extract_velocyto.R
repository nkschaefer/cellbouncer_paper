#! /usr/bin/env Rscript
library(rhdf5)

args <- commandArgs(trailingOnly=TRUE)

f <- H5Fopen(args[1])

genes <- data.frame(gid=f$row_attrs$Accession, gene=f$row_attrs$Gene)

genetab <- as.data.frame(table(genes$gene))

genes$name <- genes$gene
genes[which(genes$name %in% genetab[which(genetab$Freq > 1),]$Var1),]$name <- 
    genes[which(genes$name %in% genetab[which(genetab$Freq > 1),]$Var1),]$gid

bc <- gsub("x", "", gsub("demux:", "", f$col_attrs$CellID))

library(Matrix)

# Rows = cells
# Columns = genes
genesplice <- data.frame(gene=genes$name, spliced=colSums(f$layers$spliced))
geneunsplice <- data.frame(gene=genes$name, unspliced=colSums(f$layers$unspliced))

gene <- merge(genesplice, geneunsplice)

write.table(gene, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#cellsplice <- data.frame(barcode=bc, spliced=rowSums(f$layers$spliced))
#cellunsplice <- data.frame(barcode=bc, unspliced=rowSums(f$layers$unspliced))
#cell <- merge(cellsplice, cellunsplice)
#write.table(cell, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

