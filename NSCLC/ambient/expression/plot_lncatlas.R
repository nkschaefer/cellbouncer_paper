#! /usr/bin/env Rscript
library(pheatmap)
library(viridis)
library(reshape2)
library(dendextend)

dat <- read.table('lncAtlas_dat.txt', header=T)

dat <- dat[grep("vs", dat$Program, invert=T),]
dat$Program <- factor(dat$Program, labels=c("CellBouncer", "DecontX", "SoupX", "cellbender"),
    levels=c("cellbouncer", "decontX", "soupX", "cellbender"))

casted <- dcast(dat, Program ~ Variable, value.var="rho")

rownames(casted) <- casted$Program
casted <- casted[,-c(1)]
mat <- as.matrix(casted)

absm <- max(abs(mat))

tree_row <- hclust(dist(mat))
tree_row <- rotate(tree_row, order=c("CellBouncer", "DecontX", "SoupX", "cellbender"))

tree_col <- hclust(dist(t(mat)))
tree_col <- rotate(tree_col, order=c("Cytoplasmic_Nuclear", "Insoluble_fraction", "Membrane", "Chromatin", "Nucleolus", "Nucleoplasm"))

#tree_row <- rotate(tree_row, order=c("CellBouncer", "DecontX", "SoupX", "cellbender", "CellBouncer_vs_cellbender", "CellBouncer_vs_SoupX", "CellBouncer_vs_DecontX"))
pdf("lncAtlas_plot.pdf", width=6, height=4, bg='white')
pheatmap(mat, border_color=NA, cluster_rows=tree_row, cluster_cols=tree_col, breaks=seq(-absm, absm, (2*absm)/100), color=mako(100))
dev.off()

