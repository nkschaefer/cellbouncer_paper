#! /usr/bin/env Rscript

dx <- read.table('decontX_contam.txt', header=T)
colnames(dx)[2] <- "decontX"

wd <- read.table('demux_vcf_withdoub.contam_rate')
nd <- read.table('demux_vcf_nodoub.contam_rate')

wd <- wd[,c(1,2)]
nd <- nd[,c(1,2)]

colnames(wd) <- c("bc", "CellBouncer")
colnames(nd) <- c("bc", "CellBouncer_no_doublet")

merged <- merge(dx, wd)
merged <- merge(merged, nd)

library(reshape2)
melted <- melt(merged, id.vars=c("bc"))

library(ggplot2)

melted$variable <- factor(melted$variable, labels=c("CellBouncer_no_doublet", "CellBouncer", "decontX"),
                          levels=c("CellBouncer_no_doublet", "CellBouncer", "decontX"))

plt <- ggplot(melted) + geom_boxplot(aes(x=variable, y=value)) + 
    theme_bw() + 
    scale_y_continuous("Contamination rate") + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y=element_text(size=12),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=14))
ggsave(plt, file="ranzoni_comp.pdf", width=4, height=3.65)

