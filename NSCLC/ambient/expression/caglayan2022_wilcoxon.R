#! /usr/bin/env Rscript

prof <- read.table('contam_profs_all.txt', header=T, sep='\t')
amb <- read.table('caglayan2022_amb_genes.txt')

prof$in_set <- 0
prof[which(prof$gene %in% amb$V1),]$in_set <- 1

tab <- data.frame(prog=colnames(prof[,-c(1, length(colnames(prof)))]), U=0, p=0)

for (prog in colnames(prof[,-c(1, length(colnames(prof)))])){
    res <- wilcox.test(prof[which(prof$in_set==1),][[prog]], prof[which(prof$in_set==0),][[prog]], alternative='greater')
    tab[which(tab$prog==prog),]$U <- res$statistic
    tab[which(tab$prog==prog),]$p <- res$p.value
}

write.table(tab, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
