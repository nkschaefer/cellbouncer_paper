#! /usr/bin/env Rscript

demuxalot <- read.table('demuxalot_donor_all_500k.assignments')
demuxlet <- read.table('demuxlet_donor_all_500k.assignments')
demux_vcf <- read.table('demux_vcf_donor_all_500k.assignments')
vireo <- read.table('vireo_donor_all_500k.assignments')

cellplex <- read.table('40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_assignment_confidence_table.csv', sep=',', header=T)
cellplex$Assignment <- gsub("CMO301", "donor_1", gsub("CMO302", "donor_2", gsub("CMO303", "donor_3", 
                       gsub("CMO304", "donor_4", gsub("CMO306", "donor_5", gsub("CMO307", "donor_6",
                       gsub("CMO308", "donor_7", cellplex$Assignment)))))))
colnames(cellplex)[2] <- "donor_1"
colnames(cellplex)[3] <- "donor_2"
colnames(cellplex)[4] <- "donor_3"
colnames(cellplex)[5] <- "donor_4"
colnames(cellplex)[6] <- "donor_5"
colnames(cellplex)[7] <- "donor_6"
colnames(cellplex)[8] <- "donor_7"

cellplex$best <- apply(cellplex[,c(2,3,4,5,6,7,8)], 1, function(x){ colnames(cellplex)[which.max(x)+1] })
cellplex$second <- apply(cellplex[,c(2,3,4,5,6,7,8)], 1, function(x){ colnames(cellplex)[order(x, decreasing=T)[2] +1] })

cellplex$doub_id <- paste(cellplex$best, cellplex$second, sep="+")
cellplex[which(cellplex$second < cellplex$best),]$doub_id <- paste(cellplex[which(cellplex$second < cellplex$best),]$second,
                                                                   cellplex[which(cellplex$second < cellplex$best),]$best,
                                                                   sep='+')
cellplex$id <- cellplex$Assignment
cellplex <- cellplex[which(! cellplex$Assignment %in% c("Blanks", "Unassigned")),]
cellplex[which(cellplex$Assignment == "Multiplet"),]$id <- cellplex[which(cellplex$Assignment=="Multiplet"),]$doub_id
cellplex <- cellplex[,which(colnames(cellplex) %in% c("Barcodes", "id"))]
cellplex$Barcodes <- gsub("-1", "", cellplex$Barcodes)

colnames(cellplex) <- c("bc", "id")

demuxalot <- demuxalot[,c(1,2)]
demuxlet <- demuxlet[,c(1,2)]
demux_vcf <- demux_vcf[,c(1,2)]
vireo <- vireo[,c(1,2)]

colnames(demuxalot) <- c("bc", "id")
colnames(demuxlet) <- c("bc", "id")
colnames(demux_vcf) <- c("bc", "id")
colnames(vireo) <- c("bc", "id")

all <- rbind(cellplex, demuxalot, demuxlet, demux_vcf, vireo)
all$count <- 1
all$bc <- gsub("-1", "", all$bc)
counts <- aggregate(all$count, by=list(bc=all$bc, id=all$id), FUN=sum)
maxes <- aggregate(counts$x, by=list(bc=counts$bc), FUN=max)

colnames(counts)[3] <- "count"
colnames(maxes)[2] <- "max"

merged <- merge(counts, maxes)

maxcount <- as.data.frame(table(merged[which(merged$count==merged$max),]$bc))
ties <- unique(maxcount[which(maxcount$Freq > 1),]$Var1)

merged <- merged[which(merged$count==merged$max),]
merged <- merged[which(! merged$bc %in% ties),]

tiebreak <- cellplex[which(cellplex$bc %in% ties),]

merged <- rbind(merged[,c(1,2)], tiebreak)

write.table(merged, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)


