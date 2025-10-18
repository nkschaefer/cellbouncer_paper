#! /usr/bin/env Rscript
library(reshape2)
library(pheatmap)

ac <- read.table('mito_expr_dat.tsv', header=T)
ac <- ac[which(ac$id1 != "Chimp_CB" & ac$id2 != "Chimp_CB"),]
library(ggplot2)
anno_dat <- data.frame(id=c("Chimp", "Chimp_HC", "Human", "Human_HC", "Human_Chimp", "Chimp_Bonobo", "Bonobo_CB"),
                       genome=c("Chimp", "Human_Chimp", "Human", "Human_Chimp", "Human_Chimp", "Chimp_Bonobo", "Chimp_Bonobo"),
                       mito=c("Chimp", "Chimp", "Human", "Human", "Human_Chimp", "Chimp_Bonobo", "Bonobo"))
ord1 <- c("Chimp_Bonobo", "Human_Chimp", "Human_HC", "Chimp_HC", "Bonobo_CB", "Chimp_CB", "Human", "Chimp", "Bonobo")
ord2 <- c("Chimp+Bonobo", "Human+Chimp", "Human_HC", "Chimp_HC", "Bonobo_CB", "Chimp_CB", "Human", "Chimp", "Bonobo")
ac$id1 <- factor(ac$id1, labels=ord2, levels=ord1)
anno_dat$id <- factor(anno_dat$id, labels=ord2, levels=ord1)
anno_dat$mito <- factor(anno_dat$mito, labels=ord2, levels=ord1)
anno_dat$genome <- factor(anno_dat$genome, labels=ord2, levels=ord1)
cols1=c(Human="#7EBDC2", Chimp="#EF9943", Human_Chimp="#F3DFA2", Chimp_Bonobo="#BB4430", Bonobo="#7C3626")
colsA <- c("Human", "Chimp", "Bonobo", "Human+Chimp", "Chimp+Bonobo")
colsB <- c("#7EBDC2", "#EF9943", "#7C3626", "#2b4162", "#a63a50")
colsC <- setNames(colsB, colsA)
accpy <- merge(ac, anno_dat, by.x="id1", by.y="id")
accpy <- accpy[grep("^MT-", accpy$gene),]
accpy <- accpy[grep("MT-T", accpy$gene, invert=T),]
accpy <- accpy[grep("MT-RNR", accpy$gene, invert=T),]

plt <- ggplot(accpy) + 
    geom_vline(aes(xintercept=0), colour="lightgrey") +
    geom_vline(aes(xintercept=-1), colour="lightgrey") + 
    geom_vline(aes(xintercept=1), colour="lightgrey") +  
    geom_violin(aes(x=lfcExpr, y=id1, fill=genome, colour=mito), lwd=1.1, show.legend=FALSE) + 
    geom_text(data=anno_dat, aes(y=id, x=-11.5, label=mito, colour=mito), 
              fontface='bold', show.legend=FALSE, hjust=0, vjust=1) + 
    geom_text(data=anno_dat, aes(y=id, x=-7.5, label=genome, colour=genome), 
              fontface='bold', show.legend=FALSE, hjust=0, vjust=1) + 
    annotate("text", x=-11.5, y=8, label="MT", hjust=0, vjust=1, fontface='bold', size=12, size.unit='pt') + 
    annotate("text", x=-7.5, y=8, label="Genome", hjust=0, vjust=1, fontface='bold', size=12, size.unit='pt') + 
    scale_colour_manual(values=colsC) +
    scale_fill_manual(values=colsC) + 
    scale_x_continuous("Log2FC, Copy number + Expression", breaks=c(-2,-1,0,1,2)) + 
    theme_bw() + 
    theme(
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.line.x=element_line(linewidth=1, linetype='solid', colour='lightgrey'),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=12),
          axis.title.x=element_text(size=12))
ggsave(plt, file='lfc_expr_violin.pdf', width=5.5, height=3)

