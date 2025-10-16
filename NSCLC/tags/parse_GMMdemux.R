#! /usr/bin/env Rscript

library(stringr)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    print("Need output dir")
    q()
}

outdir <- args[1]

config <- read.table(paste(outdir, '/', 'GMM_full.config', sep=''), sep=',')
out <- read.table(paste(outdir, '/', 'GMM_full.csv', sep=''), header=T, sep=',')

colnames(config) <- c("Cluster_id", "id")
config$Cluster_id <- trimws(config$Cluster_id)
config$id <- trimws(config$id)

config$n_id <- str_count(config$id, "-")
config$type <- "S"
config[which(config$id=="negative"),]$type <- "Negative"
config[which(config$n_id==1),]$type <- "D"
config[which(config$n_id > 1),]$type <- "Multiplet"
config[which(config$type=="D"),]$id <- gsub("-", "+", config[which(config$type=="D"),]$id)

colnames(out) <- c("bc", "Cluster_id", "confidence")

merged <- merge(out, config)
merged <- merged[,c(2,4,6,3)]
merged$bc <- gsub("-1", "", merged$bc)

# Convert to 2-col
merged[which(merged$type != "S"),]$id <- "multiplet"

write.table(merged[,c(1,2)], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

