#! /usr/bin/env Rscript
library(ggplot2)
library(viridis)

profs <- read.table('contam_profs_all.txt', header=T)

rates <- read.table('contam_rates_all.txt', header=T)

dat <- data.frame(prog1=c(), prog2=c(), Correlation=c(), type=c())

for (prog1 in colnames(profs)[2:length(colnames(profs))]){
    for (prog2 in colnames(profs)[2:length(colnames(profs))]){
        if (prog1 != prog2){
            res1 <- cor(profs[[prog1]], profs[[prog2]])
            res2 <- cor(rates[[prog1]], rates[[prog2]])
            new <- data.frame(prog1=c(prog1, prog2, prog1, prog2), prog2=c(prog2, prog1, prog2, prog2), 
                 Correlation=c(res1, res1, res2, res2), type=c("genes", "genes", "cells", "cells"))
            dat <- rbind(dat, new)
        }
    }
}

levels=c("soupX", "decontX", "cellbender", "cellbouncer")
labels=c("SoupX", "DecontX", "cellbender", "CellBouncer")
dat$prog1 <- factor(dat$prog1, levels=levels, labels=labels)
dat$prog2 <- factor(dat$prog2, levels=levels, labels=labels)

print(dat)

gene <- dat[which(dat$type=="genes"),]
contam <- dat[which(dat$type=="cells"),]

gene <- gene[which(as.integer(gene$prog1) < as.integer(gene$prog2)),]
contam <- contam[which(as.integer(contam$prog1) > as.integer(contam$prog2)),]
both <- rbind(gene, contam)


#both <- read.table('contam_corr_final.tsv', header=T)

plt <- ggplot(both) + 
    geom_tile(aes(x=prog1, y=prog2, fill=Correlation)) + 
    theme_bw() + 
    geom_abline(aes(slope=1, intercept=0), lty='dotted', lwd=1.25) + 
    scale_fill_viridis(option="mako") + 
    scale_x_discrete("Cells") + 
    scale_y_discrete("Genes") + 
    theme(axis.text.y=element_text(size=12), 
          axis.title.x=element_text(size=14), 
          axis.title.y=element_text(size=14), 
          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1), 
          panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(plt, file="contam_corr.pdf", width=3.91, height=2.79)
