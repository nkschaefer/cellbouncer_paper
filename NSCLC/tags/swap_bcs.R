#! /usr/bin/env Rscript

map <- read.table('filt_barcodes_translated_map.tsv')
colnames(map) <- c("bc", "barcode")

args <- commandArgs(trailingOnly=TRUE)

counts <- read.table(args[1], header=T)

colnames(counts)[1] <- "bc"
counts <- merge(counts, map)

last <- length(colnames(counts))
counts <- counts[,c(last, seq(2, last-1))]

write.table(counts, file=args[1], sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

