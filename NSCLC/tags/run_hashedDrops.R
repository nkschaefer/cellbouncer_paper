#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("ARGS: .counts file (conf)", stderr())
    q()
}

conf <- FALSE
if (length(args) > 1 & args[2] == "conf"){
    conf <- TRUE
}

counts <- read.table(args[1], sep='\t', header=T)

suppressPackageStartupMessages(library(DropletUtils))

rownames(counts) <- counts$barcode
countsm <- t(as.matrix(counts[,-c(1)]))

res <- hashedDrops(countsm)

res$barcode <- rownames(res)

# Convert indices to names
names <- data.frame(Best=seq(1, length(rownames(countsm))), name=rownames(countsm))

res <- merge(res, names)

res[which(res$Doublet),]$name <- "multiplet"

if (conf){
    res <- res[which(res$Confident),]
}

write.table(res[,c(8,9)], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



