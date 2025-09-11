#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("ERROR: provide suffix", stderr())
    q()
}

suffix <- args[1]

species_names <- read.table(paste('demux_species_', suffix, '/SAMN19715217/species_names.txt', sep=''))
shrew_name <- paste('demux_species_', suffix, '/SAMN19715217/species_counts.txt', sep='')
rat_name <- paste('demux_species_', suffix, '/SAMN19715216/species_counts.txt', sep='')
macaque2_name <- paste('demux_species_', suffix, '/SAMN19715221/species_counts.txt', sep='')

shrew <- read.table(shrew_name)
rat <- read.table(rat_name)
#macaque1 <- read.table(macaque1_name)
macaque2 <- read.table(macaque2_name)

tot_bc <- length(rownames(shrew)) + length(rownames(rat)) + length(rownames(macaque2))

# Read randomly generated barcode list (designed to make all barcodes unique)
unique_bc <- read.table('unique_bc.txt')
unique_bc$V1 <- gsub("-1", "", unique_bc$V1)

# Add in the new (unique) barcode to each DF
start <- 1
end <- start + length(rownames(shrew))-1
shrew$bc.uniq <- unique_bc[seq(start,end),]
start <- end + 1
end <- start + length(rownames(rat))-1
rat$bc.uniq <- unique_bc[seq(start,end),]
start <- end + 1
#end <- start + length(rownames(macaque1))-1
##macaque1$bc.uniq <- unique_bc[seq(start,end),]
##start <- end + 1
end <- start + length(rownames(macaque2))-1
macaque2$bc.uniq <- unique_bc[seq(start,end),]

bc_rat_ratrhesus <- read.table('barcodes_Rat_RatRhesus.txt')
bc_rat_ratshrew <- read.table('barcodes_Rat_RatTree_Shrew.txt')
bc_rhesus_ratrhesus <- read.table('barcodes_Rhesus_RatRhesus.txt')
bc_rhesus_rhesusshrew <- read.table('barcodes_Rhesus_RhesusTree_Shrew.txt')
bc_shrew_ratshrew <- read.table('barcodes_Tree_Shrew_RatTree_Shrew.txt')
bc_shrew_rhesusshrew <- read.table('barcodes_Tree_Shrew_RhesusTree_Shrew.txt')

rat_rhesus1 <- rat[which(rat$V1 %in% bc_rat_ratrhesus$V1),]
rat_rhesus2 <- macaque2[which(macaque2$V1 %in% bc_rhesus_ratrhesus$V1),]

rat_shrew1 <- rat[which(rat$V1 %in% bc_rat_ratshrew$V1),]
rat_shrew2 <- shrew[which(shrew$V1 %in% bc_shrew_ratshrew$V1),]

rhesus_shrew1 <- macaque2[which(macaque2$V1 %in% bc_rhesus_rhesusshrew$V1),]
rhesus_shrew2 <- shrew[which(shrew$V1 %in% bc_shrew_rhesusshrew$V1),]

namevec <- c("bc", species_names$V2, 'bc_uniq')
colnames(rat_rhesus1) <- namevec
colnames(rat_rhesus2) <- namevec
colnames(rat_shrew1) <- namevec
colnames(rat_shrew2) <- namevec
colnames(rhesus_shrew1) <- namevec
colnames(rhesus_shrew2) <- namevec

# Remove these doublet cells from the singlets
rat <- rat[which(! rat$V1 %in% bc_rat_ratrhesus$V1),]
rat <- rat[which(! rat$V1 %in% bc_rat_ratshrew$V1),]
macaque2 <- macaque2[which(! macaque2$V1 %in% bc_rhesus_ratrhesus$V1),]
macaque2 <- macaque2[which(! macaque2$V1 %in% bc_rhesus_rhesusshrew$V1),]
shrew <- shrew[which(! shrew$V1 %in% bc_shrew_ratshrew$V1),]
shrew <- shrew[which(! shrew$V1 %in% bc_shrew_rhesusshrew$V1),]

# Need to scale doublets AND create doublets
# Scale Rat / Rhesus
rat_rhesus <- rat_rhesus1
scales <- rowMeans(data.frame(x=rowSums(rat_rhesus1[,-c(1,length(colnames(rat_rhesus1)))]),
               y=rowSums(rat_rhesus2[,-c(1,length(colnames(rat_rhesus2)))])))
scales1 <- rowSums(rat_rhesus1[,-c(1,length(colnames(rat_rhesus1)))]) / scales
scales2 <- rowSums(rat_rhesus2[,-c(1,length(colnames(rat_rhesus2)))]) / scales
for (coln in seq(2, length(colnames(rat_rhesus2))-1)){
    rat_rhesus1[,coln] <- round(rat_rhesus1[,coln] / scales1)
    rat_rhesus2[,coln] <- round(rat_rhesus2[,coln] / scales2)
    rat_rhesus[,coln] <- rat_rhesus1[,coln] + rat_rhesus2[,coln]
}
# Scale Rat / Shrew
rat_shrew <- rat_shrew1
scales <- rowMeans(data.frame(x=rowSums(rat_shrew1[,-c(1,length(colnames(rat_shrew1)))]),
               y=rowSums(rat_shrew2[,-c(1,length(colnames(rat_shrew2)))])))
scales1 <- rowSums(rat_shrew1[,-c(1,length(colnames(rat_shrew1)))]) / scales
scales2 <- rowSums(rat_shrew2[,-c(1,length(colnames(rat_shrew2)))]) / scales
for (coln in seq(2, length(colnames(rat_shrew1))-1)){
    rat_shrew1[,coln] <- round(rat_shrew1[,coln] / scales1)
    rat_shrew2[,coln] <- round(rat_shrew2[,coln] / scales2)
    rat_shrew[,coln] <- rat_shrew1[,coln] + rat_shrew2[,coln]
}
# Scale Rhesus / Shrew
rhesus_shrew <- rhesus_shrew1

scales <- rowMeans(data.frame(x=rowSums(rhesus_shrew1[,-c(1,length(colnames(rhesus_shrew1)))]),
                   y=rowSums(rhesus_shrew2[,-c(1,length(colnames(rhesus_shrew2)))])))
scales1 <- rowSums(rhesus_shrew1[,-c(1,length(colnames(rhesus_shrew1)))]) / scales
scales2 <- rowSums(rhesus_shrew2[,-c(1,length(colnames(rhesus_shrew2)))]) / scales
for (coln in seq(2, length(colnames(rhesus_shrew1))-1)){
    rhesus_shrew1[,coln] <- round(rhesus_shrew1[,coln] / scales1)
    rhesus_shrew2[,coln] <- round(rhesus_shrew2[,coln] / scales2)
    rhesus_shrew[,coln] <- rhesus_shrew1[,coln] + rhesus_shrew2[,coln]
}

# Store a key AND a counts table.
rat_rhesus$true <- "Rat+Rhesus"
rat_shrew$true <- "Rat+Tree_shrew"
rhesus_shrew$true <- "Rhesus+Tree_shrew"

counts_all <- rbind(rat_rhesus, rat_shrew, rhesus_shrew)
counts_all$filtered_set <- 1

bc_Rat <- read.table('barcodes_Rat.txt')
bc_Rhesus <- read.table("barcodes_Rhesus.txt")
bc_Shrew <- read.table('barcodes_Tree_Shrew.txt')

rat$true <- "Rat"
macaque2$true <- "Rhesus"
shrew$true <- "Tree_shrew"

rat$filtered_set <- 0
macaque2$filtered_set <- 0
shrew$filtered_set <- 0

rat[which(rat$V1 %in% bc_Rat$V1),]$filtered_set <- 1
macaque2[which(macaque2$V1 %in% bc_Rhesus$V1),]$filtered_set <- 1
shrew[which(shrew$V1 %in% bc_Shrew$V1),]$filtered_set <- 1

colnames(rat) <- colnames(counts_all)
colnames(macaque2) <- colnames(counts_all)
colnames(shrew) <- colnames(counts_all)

counts_all <- rbind(counts_all, rat, macaque2, shrew)

ncol <- length(colnames(counts_all))

bcmap <- counts_all[,which(colnames(counts_all) %in% c("bc", "bc_uniq"))]
key <- counts_all[,seq(ncol-2, ncol)]

counts_all <- counts_all[,c(ncol-2, seq(2,ncol-3))]

dir_out <- paste('combined_', suffix, sep='')
dir.create(dir_out, showWarnings=FALSE)

countsfile <- paste(dir_out, '/species_counts.txt', sep='')
keyfile <- paste(dir_out, '/key.txt', sep='')
bcmapfile <- paste(dir_out, '/bcmap.txt', sep='')
namesfile <- paste(dir_out, '/species_names.txt', sep='')

write.table(counts_all, file=countsfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(key, file=keyfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(bcmap, file=bcmapfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(species_names, file=namesfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

