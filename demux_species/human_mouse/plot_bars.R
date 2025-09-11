#! /usr/bin/env Rscript

bc <- read.table('../human_mouse/5k_hgmm_3p_nextgem_5k_hgmm_3p_nextgem_count_sample_filtered_barcodes.csv', sep=',')
bc$V2 <- gsub("-1", "", bc$V2)
bc$species <- ""
bc[which(bc$V1=="GRCh38"),]$species <- "Human"
bc[which(bc$V1=="GRCm39"),]$species <- "Mouse"
bc <- bc[,c(2,3)]
colnames(bc) <- c("bc", "species")

nbc <- length(rownames(bc))
nbc <- 5000
expec <- data.frame(species=c("Human", "Mouse", "Human+Mouse"), count=c(nbc*0.96*0.5, nbc*0.96*0.5, nbc*0.04))
expec$type <- "Expected"

k22 <- read.table('kmer_22_cellbouncer/species.assignments')
k22 <- k22[which(k22$V1 %in% bc$bc),]
k22 <- k22[,c(1,2)]
colnames(k22) <- c("bc", "species")
k22$count <- 1
k22agg <- aggregate(k22$count, by=list(species=k22$species), FUN=sum)
colnames(k22agg)[2] <- "count"
k22agg$type <- "kmer"

frombam <- read.table('demux_species_composite/species.assignments')
frombam <- frombam[which(frombam$V1 %in% bc$bc),]
frombam <- frombam[,c(1,2)]
frombam$V2 <- gsub("hg38", "Human", gsub("mm10", "Mouse", frombam$V2))

colnames(frombam) <- c("bc", "species")
frombam$count <- 1
frombamagg <- aggregate(frombam$count, by=list(species=frombam$species), FUN=sum)
colnames(frombamagg)[2] <- "count"
frombamagg$type <- "BAM"

bc$count <- 1
bcagg <- aggregate(bc$count, by=list(species=bc$species), FUN=sum)
colnames(bcagg)[2] <- "count"
bcagg$type <- "10X"

all <- rbind(expec, k22agg, frombamagg, bcagg)

for (t in unique(all$type)){
    all[which(all$type==t),]$count <- all[which(all$type==t),]$count / sum(all[which(all$type==t),]$count)
}

all$type <- factor(all$type, labels=c("Expected", "10X", "BAM", "k-mers"), levels=c("Expected", "10X", "BAM", "kmer"))
colnames(all)[1] <- "Species"

library(ggplot2)

plt <- ggplot(all) + geom_bar(aes(x=type, y=count, fill=Species), stat='identity') + 
    scale_y_continuous("Fraction of cells") + 
    theme_bw() + 
    scale_fill_manual(labels=c("Human", "Human+Mouse", "Mouse"), 
                      values=c("#DD6E42", "#E8DAB2", "#4F6D7A")) + 
    coord_flip() + 
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12, angle = 90, vjust = 0.5, hjust=1),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=12),
          axis.title.y=element_blank(),
          legend.position="top")

ggsave(plt, file="bars.pdf", width=6, height=3)

