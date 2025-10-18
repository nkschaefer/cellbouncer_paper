#! /usr/bin/env Rscript
library(Matrix)
library(Matrix.utils)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3){
    write("ARGS: expr_sums lib1 lib2 lib3", stderr())
    write("expr_sums is a file listing mitochondrial expression sums and whole-cell expression sums \
(see README.md)", stderr())
    write("lib1, lib2, and lib3 are the names of libraries included in the experiment", stderr())
}
expr_sum <- args[1]
libs <- read.table(args[2])
outs <- read.table(args[3], sep='\t')

spec <- read.table('bc2species.tsv', sep='\t')
colnames(spec) <- c("bc", "species")
pb_all <- NA
meta_all <- NA

mito <- read.table('mito/mito_ase_results_HC.tsv', header=T)

mito <- mito[,c(1,5)]
colnames(mito) <- c("bc", "mito")

mito2 <- spec[which(spec$species %in% c("Human_Human", "Chimp_Chimp")),]
mito2$mito <- "Human_sep"
mito2[which(mito2$species=="Chimp_Chimp"),]$mito <- "Chimp_sep"
#mito2[which(mito2$species=="Chimp_Bonobo"),]$mito <- "Chimp+Bonobo"
mito <- rbind(mito, mito2[,c(1,3)])

mito3 <- read.table('mito/mito_ase_results_CB.tsv', header=T)
mito3 <- mito3[,c(1,5)]
colnames(mito3) <- c("bc", "mito")
mito3[which(mito3$mito == "Bonobo+Chimp"),]$mito <- "Chimp_Bonobo"
mito3[which(mito3$mito == "Chimp"),]$mito <- "Chimp_CB"
mito3[which(mito3$mito == "Bonobo"),]$mito <- "Bonobo_CB"
mito <- rbind(mito, mito3)

# This is a filtered list of cells with reasonable mitochondrial expression
expr <- read.table(expr_sum, header=T)
mito <- mito[which(mito$bc %in% expr$V1),]

libs$out <- outs$V1

first <- TRUE

for (lib in args[2:length(args)]){
    write(lib, stderr())
    dirn <- libs[which(libs$V1==lib),]$out

    # Rows = genes
    # Cols = cell barcodes
    dat <- readMM(paste(lib, '_clean_mtx/matrix.mtx.gz', sep=''))
    bc <- read.table(paste(lib, '_clean_mtx/barcodes.tsv.gz', sep=""))
    gene <- read.table(paste(lib, '_clean_mtx/features.tsv.gz', sep=''))
    dat <- t(dat)
    rownames(dat) <- bc$V1
    colnames(dat) <- gene$V2
    rownames(dat) <- gsub("-1", "", rownames(dat))
    
    rownames(dat) <- paste(rownames(dat), lib, sep='-')
    dat <- dat[which(rownames(dat) %in% spec$bc),]

    meta <- data.frame(bc=rownames(dat))
    #meta <- merge(meta, spec)
    
    assn <- read.table(paste(lib '_clean.assignments', sep=''))
    assn <- assn[,c(1,2)]
    colnames(assn) <- c("bc", "id")
    meta <- merge(meta, assn)
    meta <- merge(meta, mito)

    dat <- dat[which(rownames(dat) %in% meta$bc),]
    meta <- meta[which(meta$bc %in% rownames(dat)),]
    meta <- meta[match(rownames(dat), meta$bc),]
    meta$lib <- lib

    # Aggregate
    meta$bin <- paste(meta$id, '.', lib, '.', meta$mito, sep='')
    pb <- as.data.frame(as.matrix(aggregate.Matrix(dat, meta$bin, fun='sum')))
    meta <- unique(meta[,c(5,2,3,4)])
    meta <- meta[match(rownames(pb), meta$bin),]    
    
    if (first){
        pb_all <- pb
        meta_all <- meta
        first <- FALSE
    }
    else{
        pb_all <- rbind(pb_all, pb)
        meta_all <- rbind(meta_all, meta)
    }
}

meta <- meta_all
rownames(meta) <- meta$bin

suppressPackageStartupMessages(library(DESeq2))

for (cn in colnames(pb_all)){
    pb_all[[cn]] <- as.integer(round(pb_all[[cn]]))
}

pb_all <- t(pb_all)

# Limit to just mitochondrial genes, so expression is normalized over
# the mitochondrion rather than genome-wide. This will measure post-processing
# differences rather than MT expresion/copy number differences.
pb_all <- pb_all[grep("^MT-", rownames(pb_all)),]

meta$lib_id <- paste(meta$lib, meta$id, sep='_')

dds <- DESeqDataSetFromMatrix(round(pb_all), meta, ~ lib + mito)

dds <- DESeq(dds)

saveRDS(dds, file="deseq_results_mito_norm.rds")

