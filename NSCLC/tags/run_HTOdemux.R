#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("ARGS: tags_mtx", stderr())
    q()
}

suppressPackageStartupMessages(library(Seurat))

tags <- ReadMtx(paste(args[1], '/matrix.mtx.gz', sep=''), 
    paste(args[1], '/barcodes.tsv.gz', sep=''), 
    paste(args[1], '/features.tsv.gz', sep=''),
    feature.column=1)

tags <- tags + 1

bc <- read.table(paste(args[1], '/barcodes.tsv.gz', sep=''))

tags.seurat <- CreateSeuratObject(counts=tags, sparse=T)

# Need to trick Seurat into thinking there's gene expression data too
m_empty <- t(matrix(NA, length(rownames(bc)), 0, dimnames=list(bc$V1, c())))

gex.seurat <- CreateSeuratObject(counts=as(m_empty, "CsparseMatrix"))

gex.seurat[["HTO"]] <- CreateAssayObject(counts=tags)

# Normalize
gex.seurat <- NormalizeData(gex.seurat, assay="HTO", normalization.method="CLR")

gex.seurat <- HTODemux(gex.seurat, assay="HTO", positive.quantile=0.99)

df <- as.data.frame(gex.seurat$HTO_classification)
df$barcode <- rownames(df)
df <- df[,c(2,1)]
df$type <- gex.seurat$HTO_classification.global

colnames(df) <- c("barcode", "id", "type")

df[which(df$type != "Singlet"),]$id <- "multiplet"

write.table(df[,c(1,2)], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

