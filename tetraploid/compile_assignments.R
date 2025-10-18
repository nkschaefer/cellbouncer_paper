#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
    write("ARGS: demux_vcf.assignments demux_tags.assignments mito_expr_genes", stderr())
    write("Both files should be concatenated across all libraries, with library IDs appended to barcodes",
          stderr())
    q()
}

bcfilt <- read.table(args[3], sep='\t')

vcf <- read.table(args[1], sep='\t')
colnames(vcf) <- c("bc", "vcf", "type", "llr")

ms <- read.table(args[2], sep='\t')
ms <- ms[which(ms$V3=="S"),]
ms <- ms[,c(1,2)]
colnames(ms) <- c("bc", "ms")

merged <- merge(vcf, ms)
merged <- merged[which(merged$bc %in% bcfilt$V1),]

msmap <- read.table('multiseq/msid2cellines.tsv')
colnames(msmap) <- c("ms", "cl")

merged <- merge(merged, msmap)
merged <- merged[which(merged$vcf == merged$cl),]

# Load ID -> species mapping
ids <- read.table('id2species.tsv', sep='\t')
colnames(ids) <- c("vcf", "species")

merged <- merge(merged, ids)

b2s <- merged[,which(colnames(merged) %in% c("bc", "species"))]
write.table(b2s, file="bc2species.tsv", sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

b2i <- merged[,which(colnames(merged) %in% c("bc", "vcf"))]
write.table(b2i, file="bc2id.tsv", sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# Save data separated by library, and library+species.
merged$lib <- apply(merged, 1, function(x){ strsplit(x[1], '-')[[1]][2] })
for (lib in unique(merged$lib)){
    subs <- merged[which(merged$lib==lib),which(colnames(merged) %in% c("bc", "vcf", "type", "llr"))]
    write.table(subs, file=paste(lib, "_clean.assignments", sep=''), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

    hh <- merged[which(merged$lib == lib & merged$species=="Human_Human"),which(colnames(merged) %in% c("bc", "vcf", "type", "llr"))]
    cc <- merged[which(merged$lib == lib & merged$species=="Chimp_Chimp"),which(colnames(merged) %in% c("bc", "vcf", "type", "llr"))]
    hc <- merged[which(merged$lib == lib & merged$species=="Human_Chimp"),which(colnames(merged) %in% c("bc", "vcf", "type", "llr"))]
    cb <- merged[which(merged$lib == lib & merged$species=="Chimp_Bonobo"),which(colnames(merged) %in% c("bc", "vcf", "type", "llr"))]

    write.table(hh, file=paste(lib, "_hh.assignments", sep=''), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(cc, file=paste(lib, "_cc.assignments", sep=''), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(hc, file=paste(lib, "_hc.assignments", sep=''), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(cb, file=paste(lib, "_cb.assignments", sep=''), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
}
