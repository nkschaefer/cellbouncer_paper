#! /usr/bin/env Rscript
library(ggplot2)
library(viridis)

dat <- read.table('contam_profs_all.txt', header=T)
splice <- read.table('gene_splice.txt', header=T)
ll <- read.table('bulkprops_ll.txt')
ll <- ll[,c(2,3)]
colnames(ll) <- c("gene", "ll")

llF <- ecdf(ll$ll)
ll$ll.p <- llF(ll$ll)

merged <- merge(dat, ll)
merged <- merge(merged, splice)

plt <- ggplot(merged) + 
    geom_point(aes(x=cellbouncer, y=ll.p, colour=spliced+1)) + 
    scale_x_log10("Fraction of ambient RNA") + 
    scale_y_continuous("Genotype log likelihood percentile") + 
    scale_colour_viridis(name="Spliced", option="mako", trans='log10') + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=14),
          legend.title=element_text(size=14),
          legend.text=element_text(size=12))

ggsave(plt, file="cb_vs_ll_spliced.pdf", width=6, height=3.5)

