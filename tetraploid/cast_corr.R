#! /usr/bin/env Rscript

corrmap <- read.table('flank_corr.tsv')
library(reshape2)

colnames(corrmap)[1] <- "gene1"
colnames(corrmap)[2] <- "gene2"

casted <- dcast(corrmap, gene1 + gene2 ~ V7, value.var="V4")

write.table(casted, file="flank_corr_casted.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

