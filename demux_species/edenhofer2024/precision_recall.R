#! /usr/bin/env Rscript

a1 <- read.table('demux_species_out/combined/species.assignments')
a2 <- read.table('edenhofer2024_demux_species/species.assignments')

colnames(a1) <- c("bc", "id", "type", "llr")
colnames(a2) <- c("bc", "id2", "type2", "llr2")

merged <- merge(a1, a2)

prec.H <- length(rownames(merged[which(merged$id=="Human" & (merged$id2 == "Human" | merged$id2 == "Both" | merged$id2 == "Both+Human")),]))
prec.M <- length(rownames(merged[which(merged$id=="Macaque" & (merged$id2 == "Macaque" | merged$id2 == "Both" | merged$id2 == "Both+Macaque")),]))
prec.B <- length(rownames(merged[which(merged$id=="Human+Macaque" & (merged$id2 == "Both" | merged$id2 == "Both+Human" | merged$id2 == "Both+Macaque" | merged$id2 == "Human+Macaque")),]))
precision <- (prec.H + prec.M + prec.B)/length(rownames(merged))

recall.H <- length(rownames(merged[which(merged$id2=="Human" & merged$id=="Human"),]))
recall.M <- length(rownames(merged[which(merged$id2=="Macaque" & merged$id=="Macaque"),]))
recall <- (recall.H + recall.M)/(length(rownames(merged[which(merged$id2 %in% c("Human", "Macaque")),])))

f1 <- (2*precision*recall)/(precision + recall)

df <- data.frame(precision=precision, recall=recall, f1=f1)
write.table(df, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

counts <- read.table('edenhofer2024_demux_species/species_counts.txt')
colnames(counts) <- c("bc", "Human", "Macaque", "Both")

rownames(counts) <- counts$bc
counts <- counts[,-c(1)]
counts <- as.matrix(counts)
counts <- counts[which(rowSums(counts) > 10),]
counts <- counts / rowSums(counts)

library(pheatmap)
library(viridis)

meta <- data.frame(row.names=a1$bc, species=a1$id)

cols <- c("#FE4A49", "#2AB7CA", "#FED766")
cols <- setNames(cols, c("Human", "Macaque", "Human+Macaque"))
cols <- list(species=cols)

pdf("edenhofer_species_heatmap.pdf", bg='white', height=3.5, width=6)
pheatmap(counts, 
         cluster_cols=FALSE, 
         show_rownames=FALSE, 
         annotation_row=meta, 
         annotation_colors=cols,
         breaks=seq(0,1, 0.01), 
         color=viridis(100, option='mako'))
dev.off()

