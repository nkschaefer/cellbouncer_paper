#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("Please provide path to filtered_feature_bc_matrix", stderr())
    q()
}
filt <- args[1]

lnc <- read.table('lncAtlas.txt', header=T)
lnc <- lnc[which(!is.na(lnc$Value)),]

lncmean <- aggregate(lnc$Value, by=list(gid=lnc$gid, location=lnc$location), FUN=mean)
colnames(lncmean)[3] <- "RCI"

profs <- read.table('contam_profs_all.txt', header=T)

featurefile <- paste(filt, '/features.tsv.gz', sep='')
features <- read.table(featurefile, sep='\t')
features <- features[which(features$V2 %in% profs$gene),]

tab <- as.data.frame(table(features$V2))
features <- features[which(features$V3=="Gene Expression"),]
features <- features[which(! features$V2 %in% tab[which(tab$Freq > 1),]$Var1),]
features <- features[,c(1,2)]
colnames(features) <- c("gid", "gene")

profs <- merge(profs, features)
profs <- merge(profs, lncmean)

locs <- unique(profs$location)
progs <- c("cellbouncer", "soupX", "decontX", "cellbender")
other_progs <- c("soupX", "decontX", "cellbender")

corr <- data.frame(Program=c(), Variable=c(), rho=c())

for (loc in locs){
    for (prog in progs){
        res <- cor(profs[which(profs$location==loc),][[prog]],
                   profs[which(profs$location==loc),]$RCI,
                   method='spearman',
                   use='pairwise.complete.obs')
        corr <- rbind(corr, data.frame(Program=c(prog), Variable=c(loc), rho=c(res)))
    }
    for (other_prog in other_progs){
        res <- cor(profs[which(profs$location==loc),]$cellbouncer / (
                     profs[which(profs$location==loc),]$cellbouncer + profs[which(profs$location==loc),][[other_prog]]),
                   profs[which(profs$location==loc),]$RCI,
                   method='spearman',
                   use='pairwise.complete.obs')
        corr <- rbind(corr, data.frame(Program=c(paste("cellbouncer_vs_", other_prog, sep='')),
                                        Variable=c(loc), rho=c(res)))
    }
}

write.table(corr, file='lncAtlas_dat.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

