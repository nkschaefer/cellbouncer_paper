#! /usr/bin/env Rscript

library(demuxmix)

args <- commandArgs(trailingOnly=TRUE)
counts <- read.table(args[1], header=T)

rownames(counts) <- counts$barcode
counts <- counts[,-c(1)]

# This is stupid that it can't handle this, and we should complain about it in the paper
counts <- counts[,which(! colnames(counts) %in% c("CMO305", "CMO309", "CMO310", "CMO311", "CMO312"))]

counts <- t(as.matrix(counts))

sink(nullfile())
res <- demuxmix(counts, model="naive")
dat <- as.data.frame(dmmClassify(res))
sink()

dat$barcode <- rownames(dat)
dat <- dat[which(dat$Type != "uncertain"),]
dat <- dat[which(dat$Type != "negative"),]
dat[which(dat$Type=="singlet"),]$Type <- "S"
dat[which(dat$Type=="multiplet"),]$HTO <- "multiplet"
dat[which(dat$Type=="multiplet"),]$Type <- "D"

dat <- dat[,c(4,1,3,2)]

dat <- dat[,c(1,2)]

write.table(dat, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

