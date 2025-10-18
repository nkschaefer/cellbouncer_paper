#! /usr/bin/env Rscript
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2){
    write("ARGS: mito_expr_gene mito_expr_gene_sum", stderr())
    q()
}
expr <- read.table(args[1], header=T)

spec <- read.table('bc2species.tsv')
colnames(spec) <- c("bc", "species")
casted <- dcast(expr, bc ~ gene, value.var="count")
casted <- merge(casted, spec)

exprgene <- expr

id <- read.table('mito/mito_ase_results_HC.tsv', header=T)
id2 <- read.table('mito/mito_ase_results_CB.tsv', header=T)
id2[which(id2$id=="Chimp"),]$id <- "Chimp_CB"
id2[which(id2$id=="Bonobo"),]$id <- "Bonobo_CB"

id <- id[,c(1,5)]
id2 <- id2[,c(1,5)]
id <- rbind(id, id2)

casted <- merge(casted, id, all.x=TRUE)
casted[which(casted$species=="Human_Human"),]$id <- "Human"
casted[which(casted$species=="Chimp_Chimp"),]$id <- "Chimp"
casted <- casted[which(!is.na(casted$id)),]

expr <- read.table(args[2], header=T)

casted <- casted[which(casted$bc %in% expr[which(expr$mito/expr$tot < 0.2),]$bc),]

gp <- data.frame(gene1=c(), gene2=c(), corr=c(), key=c())
genes <- colnames(casted)[which(! colnames(casted) %in% c("bc", "species", "id"))]
casted$key <- ""
casted[which(casted$species=="Human_Human"),]$key <- "Human"
casted[which(casted$species=="Chimp_Chimp"),]$key <- "Chimp"
casted[which(casted$species=="Human_Chimp" & casted$id == "Human"),]$key <- "Human_HC"
casted[which(casted$species=="Human_Chimp" & casted$id == "Chimp"),]$key <- "Chimp_HC"
casted[which(casted$species=="Human_Chimp" & casted$id == "Chimp+Human"),]$key <- "Human_Chimp"
casted[which(casted$id=="Bonobo+Chimp"),]$key <- "Chimp_Bonobo"
casted[which(casted$id=="Bonobo_CB"),]$key <- "Bonobo_CB"
casted[which(casted$id=="Chimp_CB"),]$key <- "Chimp_CB"

for (key in unique(casted$key)){
    sub <- casted[which(casted$key==key),]
    for (i1 in seq(1, length(genes)-1)){
        g1 <- genes[i1]
        for (i2 in seq(i1+1, length(genes))){
            g2 <- genes[i2]
            row <- data.frame(gene1=g1, gene2=g2, corr=cor(sub[[g1]], sub[[g2]], use='pairwise.complete.obs'), key=key)
            gp <- rbind(gp, row)
        }
    }
}


gpb <- gp[,c(2,1,3,4)]
colnames(gpb)[1] <- "gene1"
colnames(gpb)[2] <- "gene2"
gp <- rbind(gp, gpb)
gtf <- read.table('mito/mito.bed')
gp$gene1 <- factor(gp$gene1, labels=gtf$V4, levels=gtf$V4)
gp$gene2 <- factor(gp$gene2, labels=gtf$V4, levels=gtf$V4)

gpc <- dcast(gp, gene1 + gene2 ~ key, value.var='corr')
gpc <- gpc[grep("MT-T", gpc$gene1, invert=T),]
gpc <- gpc[grep("MT-T", gpc$gene2, invert=T),]
gpc <- gpc[grep("MT-RNR", gpc$gene1, invert=T),]
gpc <- gpc[grep("MT-RNR", gpc$gene2, invert=T),]

write.table(gpc, file="flank_corr_dat.tsv", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

