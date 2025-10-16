#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

perc <- as.numeric(args[1])

if (perc <= 0 | perc >= 1){
    write("ERROR: perc must be between 0 and 1", stderr())
    q()
}
counts <- read.table('full.counts', header=T)
true <- read.table('../ids_true_consensus.txt')
colnames(true) <- c("barcode", "id")
true$type <- "S"
true[grep("+", true$id, fixed=T),]$type <- "D"
true$id1 <- apply(true, 1, function(x){ strsplit(x[2], '+', fixed=T)[[1]][1] })
true$id2 <- apply(true, 1, function(x){ strsplit(x[2], "+", fixed=T)[[1]][2] })

merged <- merge(counts, true)

library(reshape2)
melted <- melt(merged, id.vars=c("barcode", "id", "type", "id1", "id2"))

perc <- 1 - perc
melted[which(melted$id1 == melted$variable),]$value <- round(perc*melted[which(melted$id1 == melted$variable),]$value)
melted[which(melted$id2 == melted$variable),]$value <- round(perc*melted[which(melted$id2 == melted$variable),]$value)

casted <- dcast(melted[,c(1,6,7)], barcode ~ variable, value.var='value')
write.table(casted, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

