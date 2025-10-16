#! /usr/bin/env Rscript
library(cellhashR)

args <- commandArgs(trailingOnly=TRUE)

counts <- read.table(args[1], header=T)
rownames(counts) <- counts$barcode
counts <- t(as.matrix(counts[,-c(1)]))

# This is another program that prints stupid stuff to stdout, which is bad and dumb.
sink(nullfile())
# THIS IS SLOW AND WILL FORCE YOU TO LOOK AT A HUNDRED PLOTS YOU DON'T CARE ABOUT
x <- GenerateCellHashingCalls(counts, methods=c("bff_cluster"), doTSNE=FALSE, doHeatmap=FALSE)
sink()

df <- as.data.frame(x)

df$bff_cluster <- as.character(df$bff_cluster)
df <- df[which(df$consensuscall.global != "Negative"),]
df[which(df$consensuscall.global == "Doublet"),]$bff_cluster <- "multiplet"

write.table(df[,c(1,2)], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

