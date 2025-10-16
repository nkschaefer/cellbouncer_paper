#! /usr/bin/env Rscript

true <- read.table('../ids_true_consensus.txt', sep='\t')

true[grep("+", true$V2, fixed=T),]$V2 <- "multiplet"

f <- file('stdin')

test <- read.table(f, sep="\t")

colnames(test) <- c("barcode", "id")
colnames(true) <- c("barcode", "true")

merged <- merge(true, test, all.x=TRUE, all.y=FALSE)

match <- length(rownames(merged[which(merged$id==merged$true),]))
prec <- match/length(rownames(merged[which(!is.na(merged$id)),]))
recall <- match / length(rownames(merged[which(!is.na(merged$true)),]))

f <- 2*(prec*recall)/(prec + recall)

prec_S <- length(rownames(merged[which(merged$id==merged$true & merged$id != "multiplet"),])) / 
    length(rownames(merged[which(merged$id != "multiplet"),]))
prec_D <- length(rownames(merged[which(merged$id==merged$true & merged$id == "multiplet"),])) / 
    length(rownames(merged[which(merged$id == "multiplet"),]))

recall_S <- length(rownames(merged[which(merged$id==merged$true & merged$true != "multiplet"),])) / 
    length(rownames(merged[which(merged$true != "multiplet"),]))
recall_D <- length(rownames(merged[which(merged$id== merged$true & merged$true == "multiplet"),])) / 
    length(rownames(merged[which(merged$true == "multiplet"),]))

f_S <- 2*(prec_S*recall_S)/(prec_S + recall_S)
f_D <- 2*(prec_D*recall_D)/(prec_D + recall_D)

write(paste(prec, recall, f, prec_S, recall_S, f_S, prec_D, recall_D, f_D, sep='\t'), stdout())

