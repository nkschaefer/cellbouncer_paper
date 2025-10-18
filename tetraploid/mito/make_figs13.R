#! /usr/bin/env Rscript
library(ggplot2)

haps <- read.table('mito.haps2')
ids <- read.table('mito.ids')
rownames(haps) <- ids$V1

d <- as.matrix(dist(haps))
d <- as.data.frame(d)
d$hap1 <- rownames(d)
library(reshape2)
d <- melt(d, id.vars=c("hap1"))
colnames(d) <- c("hap1", "hap2", "dist")

within <- read.table('intra_species_dat.tsv', sep='\t', header=T)
within <- within[which(within$species %in% c("Human_Human", "Chimp_Chimp")),]

s1 <- aggregate(within$count1, by=list(bc=within$bc, hap1=within$id1, hap2=within$id2, species=within$species), FUN=sum)
s2 <- aggregate(within$count2, by=list(bc=within$bc, hap1=within$id1, hap2=within$id2, species=within$species), FUN=sum)

colnames(s1)[5] <- "count1"
colnames(s2)[5] <- "count2"
merged <- merge(s1, s2)

agg <- aggregate(abs(0.5 - merged$count1/(merged$count1+merged$count2)),
                 by=list(hap1=merged$hap1, hap2=merged$hap2, species=merged$species), 
                 FUN=mean)
colnames(agg)[4] <- "ratio"
merged <- merge(agg, d)

plt <- ggplot(merged) + 
    geom_point(aes(x=dist, y=ratio, colour=species), size=2) + 
    theme_bw() + 
    scale_x_continuous("Euclidean distance between MT haplotypes") + 
    scale_y_continuous("Mean abs(difference from 0.5)") + 
    ggtitle("MT allele fraction: within-species composite cells") + 
    scale_colour_discrete("Species") + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))
ggsave(plt, file="figs13.pdf", width=5, height=3.5)
ggsave(plt, file="figs13.png", width=5, height=3.5)


