#! /usr/bin/env Rscript

hc <- read.table('hc_data.txt', sep='\t')
cb <- read.table('cb_data.txt', sep='\t')
hh <- read.table('hh_data.txt', sep='\t')
cc <- read.table('cc_data.txt', sep='\t')

colnames(hc) <- c("bc", "s1", "s2", "gene", "count1", "count2", "lib")
colnames(cb) <- c("bc", "s1", "s2", "gene", "count1", "count2", "lib")

colnames(hh) <- c("bc", "id1", "id2", "gene", "count1", "count2", "lib")
colnames(cc) <- c("bc", "id1", "id2", "gene", "count1", "count2", "lib")

hh$species <- "Human_Human"
cc$species <- "Chimp_Chimp"
hc$species <- "Human_Chimp"
cb$species <- "Chimp_Bonobo"

inter <- rbind(hc, cb)
intra <- rbind(hh, cc)

# Assume lib names have already been appended to barcodes
#inter$bc <- paste(inter$bc, inter$lib, sep='-')
#intra$bc <- paste(intra$bc, intra$lib, sep='-')

write.table(inter, file="inter_species_dat.tsv", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(intra, file="intra_species_dat.tsv", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)


