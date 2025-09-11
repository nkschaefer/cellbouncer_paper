#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("Provide output directory", stderr())
    q()
}

dirn <- args[1]

d1 <- read.table(paste(dirn, 'lane1_SI_TT_G1/species_counts.txt', sep='/'))
d2 <- read.table(paste(dirn, 'lane1_SI_TT_H1/species_counts.txt', sep='/'))

colnames(d1) <- c("bc", "x1", "y1")
colnames(d2) <- c("bc", "x2", "y2")

merged <- merge(d1, d2, all.x=TRUE, all.y=TRUE)

merged[which(is.na(merged$x1)),]$x1 <- 0
merged[which(is.na(merged$x2)),]$x2 <- 0
merged[which(is.na(merged$y1)),]$y1 <- 0
merged[which(is.na(merged$y2)),]$y2 <- 0
merged$x <- merged$x1 + merged$x2
merged$y <- merged$y1 + merged$y2

merged <- merged[,c(1,6,7)]

write.table(merged, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

