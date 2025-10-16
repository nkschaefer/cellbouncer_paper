#! /usr/bin/env Rscript

decontX_contam <- read.table('decontX_contam.txt', header=T, sep='\t')
decontX_genes <- read.table('decontX_genes.txt', sep='\t')
soupX_genes <- read.table('soupX_genes.txt', sep='\t')
soupX_contam <- read.table('soupX_contam.txt', sep='\t')
cellbender_contam <- read.table('cellbender_contam.txt', sep='\t')
cellbender_genes <- read.table('cellbender_genes.txt', sep='\t')
cellbouncer_contam <- read.table('demux_vcf_donor_all_500k.contam_rate')
cellbouncer_genes <- read.table('demux_vcf_donor_all_500k.gex_profile', header=T, sep='\t')

colnames(decontX_genes) <- c("gene", "decontX")
colnames(decontX_contam) <- c("bc", "decontX")

soupX_genes <- soupX_genes[,c(1,2)]
colnames(soupX_genes) <- c("gene", "soupX")
soupX_contam <- soupX_contam[,c(1,5)]
colnames(soupX_contam) <- c("bc", "soupX")

colnames(cellbender_contam) <- c("bc", "cellbender")
colnames(cellbender_genes) <- c("gene", "cellbender")

cellbouncer_contam <- cellbouncer_contam[,c(1,2)]
colnames(cellbouncer_contam) <- c("bc", "cellbouncer")
cellbouncer_genes <- cellbouncer_genes[,c(1,2)]
colnames(cellbouncer_genes) <- c("gene", "cellbouncer")

rates <- merge(decontX_contam, soupX_contam)
rates <- merge(rates, cellbender_contam)
rates <- merge(rates, cellbouncer_contam)

profs <- merge(decontX_genes, soupX_genes)
profs <- merge(profs, cellbender_genes)
profs <- merge(profs, cellbouncer_genes)

write.table(rates, file="contam_rates_all.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(profs, file="contam_profs_all.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

