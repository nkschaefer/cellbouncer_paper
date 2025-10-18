#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))
all <- readRDS('deseq_results.rds')
mt <- readRDS('deseq_results_mito_norm.rds')

mitos <- c("Chimp", "Human", "Chimp_sep", "Human_sep", "Chimp.Human", "Chimp_Bonobo", "Chimp_CB", "Bonobo_CB")
dats <- data.frame()
first <- TRUE

get_dats <- function(dat, kw){
    dats <- data.frame()
    first <- TRUE
    for (i1 in seq(1, length(mitos)-1)){
        m1 <- mitos[i1]
        for (i2 in seq(i1+1, length(mitos))){
            m2 <- mitos[i2]
            res <- results(dat, c("mito", m1, m2))
            m1b <- gsub("_sep", "sep", m1)
            m2b <- gsub("_sep", "sep", m2)
            m1b <- gsub("Chimp_Bonobo", "ChimpBonobo", m1b)
            m2b <- gsub("Chimp_Bonobo", "ChimpBonobo", m2b)
            m1b <- gsub("Chimp_CB", "ChimpCB", m1b)
            m2b <- gsub("Chimp_CB", "ChimpCB", m2b)
            m1b <- gsub("Bonobo_CB", "BonoboCB", m1b)
            m2b <- gsub("Bonobo_CB", "BonoboCB", m2b)
            res <- res[,c(2,3,6)]
            colnames(res)[1] <- paste("lfc", m1b, m2b, kw, sep='_')
            colnames(res)[2] <- paste("se", m1b, m2b, kw, sep='_')
            colnames(res)[3] <- paste('p', m1b, m2b, kw, sep='_')
            res <- as.data.frame(res)
            if (first){
                dats <- res
                first <- FALSE
            } else{
                dats <- cbind(dats, res)
            }
        }
    }
    return(dats)
}

d1 <- get_dats(all, 'all')
d2 <- get_dats(mt, 'mt')

d1$gene <- rownames(d1)
d2$gene <- rownames(d2)

library(reshape2)
melt1 <- melt(d1, id.vars=c("gene"))
melt2 <- melt(d2, id.vars=c("gene"))
melted <- rbind(melt1, melt2)
melted$metric <- apply(melted, 1, function(x){ strsplit(x[2], '_')[[1]][1] })
melted$comp <- apply(melted, 1, function(x){ splits <- strsplit(x[2], '_')[[1]]
    paste(splits[2], splits[3], sep='_') })
melted$norm <- apply(melted, 1, function(x){ strsplit(x[2], '_')[[1]][4] })

all <- dcast(melted[which(melted$norm=="all"),], gene + comp ~ metric, value.var="value")
mt <- dcast(melted[which(melted$norm=="mt"),], gene + comp ~ metric, value.var='value')

merged <- merge(all, mt, by=c("gene", "comp"), all.x=TRUE)
merged$lfcDiff <- merged$lfc.x - merged$lfc.y
merged$seDiff <- sqrt(merged$se.x^2 + merged$se.y^2)
merged$waldDiff <- merged$lfcDiff^2 / merged$seDiff^2
merged$pDiff <- pchisq(merged$waldDiff, 1, lower.tail=FALSE)
merged$pDiffAdj <- 0
for (comp in unique(merged$comp)){
    merged[which(merged$comp==comp),]$pDiffAdj <- p.adjust(merged[which(merged$comp==comp),]$pDiff, method='fdr')
}
merged$id1 <- apply(merged, 1, function(x){ strsplit(x[2], '_')[[1]][1] })
merged$id2 <- apply(merged, 1, function(x){ strsplit(x[2], '_')[[1]][2] })
res <- data.frame(gene=merged$gene, id1=merged$id1, id2=merged$id2, 
                  lfc=merged$lfc.x, se=merged$se.x, p=merged$p.x,
                  lfcExpr=merged$lfcDiff, pExpr=merged$pDiffAdj,
                  seExpr=merged$seDiff, lfcProc=merged$lfc.y, 
                  pProc=merged$p.y, seProc=merged$se.y)
res2 <- res
res2$lfcExpr <- -res2$lfcExpr
res2$lfcProc <- -res2$lfcProc
res2$lfc <- -res2$lfc
res2 <- res2[,c(1,3,2,4,5,6,7,8,9,10,11,12)]
colnames(res2)[2] <- "id1"
colnames(res2)[3] <- "id2"
res <- rbind(res, res2)
res[which(res$id1=="Chimp"),]$id1 <- "Chimp_HC"
res[which(res$id2=="Chimp"),]$id2 <- "Chimp_HC"
res[which(res$id1=="Chimpsep"),]$id1 <- "Chimp"
res[which(res$id2=="Chimpsep"),]$id2 <- "Chimp"
res[which(res$id1=="Human"),]$id1 <- "Human_HC"
res[which(res$id2=="Human"),]$id2 <- "Human_HC"
res[which(res$id1=="Humansep"),]$id1 <- "Human"
res[which(res$id2=="Humansep"),]$id2 <- "Human"
res[which(res$id1=="Chimp.Human"),]$id1 <- "Human_Chimp"
res[which(res$id2=="Chimp.Human"),]$id2 <- "Human_Chimp"
res[which(res$id1=="ChimpBonobo"),]$id1 <- "Chimp_Bonobo"
res[which(res$id2=="ChimpBonobo"),]$id2 <- "Chimp_Bonobo"
res[which(res$id1=="ChimpCB"),]$id1 <- "Chimp_CB"
res[which(res$id2=="ChimpCB"),]$id2 <- "Chimp_CB"
res[which(res$id1=="BonoboCB"),]$id1 <- "Bonobo_CB"
res[which(res$id2=="BonoboCB"),]$id2 <- "Bonobo_CB"

write.table(res, file="mito_expr_dat.tsv", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

