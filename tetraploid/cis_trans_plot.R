#! /usr/bin/env Rscript
library(ggplot2)
library(ggthemes)
library(ggsci)

dat <- read.table('flank_corr_casted.tsv', header=T)

dat$Net.Change <- log2(dat$Human) - log2(dat$Chimp)
dat$MT <- log2(dat$Human_HC) - log2(dat$Chimp_HC)
dat$transH <- (log2(dat$Chimp_HC) - log2(dat$Chimp))
dat$transC <- (log2(dat$Human_HC) - log2(dat$Human))
dat$nuclear <- dat$transH - dat$transC

library(reshape2)
melted <- melt(dat[,which(colnames(dat) %in% c("gene1", "gene2", "Net.Change", "MT", "nuclear"))], 
               id.vars=c("gene1", "gene2"))

melted$type <- "Mutation type"
melted[which(melted$variable=="Net.Change"),]$type <- "Net Change"

cols=c(MT="#817E9F", nuclear="#7FC29B", Net.Change="#F3DFA2")

melted$gene_pair <- paste(melted$gene1, melted$gene2)
gpord <- c("MT-ND1 MT-ND2", "MT-ND2 MT-CO1", "MT-CO1 MT-CO2",
           "MT-CO2 MT-ATP8", "MT-ATP8 MT-ATP6", "MT-ATP6 MT-CO3",
           "MT-CO3 MT-ND3", "MT-ND3 MT-ND4L", "MT-ND4L MT-ND4",
           "MT-ND4 MT-ND5", "MT-ND5 MT-CYB")
gpord <- rev(gpord)
melted$gene_pair <- factor(melted$gene_pair, labels=gpord, levels=gpord)

mv <- max(c(abs(dat$MT + dat$nuclear), abs(dat$MT), abs(dat$nuclear)))

plt <- ggplot(melted) + 
    geom_bar(aes(x=value, fill=variable, y=gene_pair), stat='identity') + 
    facet_grid(.~type) + 
    scale_fill_manual(values=cols) + 
    theme_bw() + 
    scale_x_continuous("Log2(Human corr / Chimp corr)", limits=c(-mv,mv)) + 
    scale_y_discrete("Gene pair") + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.text.x=element_text(size=12, angle=90, vjust=0.5, hjust=1),
          axis.title.y=element_text(size=14),
          strip.background=element_blank(),
          strip.text.x=element_text(size=12, face='bold'))
ggsave(plt, file="cis_trans_proc.pdf", width=6.5, height=3.5)

