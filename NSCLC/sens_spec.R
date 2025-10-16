#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

prog <- ""
donor <- ""
num <- ""

split <- strsplit(args[1], '_')[[1]]
if (split[1] == "demux" & split[2] == "vcf"){
    prog <- "demux_vcf"
    if (split[3] == "donor"){
        donor <- paste(split[3], split[4], sep='_')
        num <- split[5]
    } else{
        donor <- split[3]
        num <- split[4]
    }
} else if (split[1] == "demux" & split[2] == "vcf2"){
    prog <- "demux_vcf2"
    if (split[3] == "donor"){
        donor <- paste(split[3], split[4], sep="_")
        num <- split[5]
    } else{
        donor <- split[3]
        num <- split[4]
    }   
}else{
    prog <- split[1]
    if (split[2] == "donor"){
        donor <- paste(split[2], split[3], sep="_")
        num <- split[4]
    } else{
        donor <- split[2]
        num <- split[3]
    }
}
num <- gsub(".assignments", "", num)
num <- gsub("k", "", num)

if ( prog == "demux_vcf") {
    prog <- "CellBouncer"
} else if (prog == "demux_vcf2"){
    prog <- "CellBouncer_contam"
}
test <- read.table(args[1], sep="\t")
test$V1 <- gsub("-1", "", test$V1)

colnames(test) <- c("bc", "id", "type", "llr")

#true <- read.table('chrM_nofilter_withdoublet.assignments2')
#true <- true[,c(1,2)]

true <- read.table('ids_true_consensus.txt')
#true <- read.table('ids_true_notie.txt')
#true <- read.table("ids_true.txt")
#true <- true[,c(1,3)]

colnames(true) <- c("bc", "id.true")
true$type.true <- "S"
true[grep("+", true$id.true, fixed=T),]$type.true <- "D"
if (donor == "multi"){
    true <- true[which(true$type.true=="D"),]
} else if (donor == "donor_all"){
    # Do not subset  
    x <- 1
} else if (length(grep("ds", donor)) > 0){
    # Do not subset
    x <- 1
} else if (donor == "sc"){
    # leave alone 
    n <- gsub(".assignments", "", split[3])
    if (prog == "CellBouncer"){
        n <- gsub(".assignments", "", split[4])
    } else if (prog == "CellBouncer_contam"){
        n <- gsub(".decontam.assignments", "", split[4])
        num <- gsub(".decontam", "", num)
    }
    n <- gsub(".decontam", "", n)
    bcs <- read.table(paste('samp', n, '.barcodes', sep=''))
    true <- true[which(true$bc %in% bcs$V1),]
}else{
    
    true <- true[which(true$id.true==donor),]
}

num <- gsub(".decontam", "", num)

merged <- merge(test, true, all.y=TRUE)

meanllr <- mean(merged$llr, na.rm=TRUE)
prec_all <- length(rownames(merged[which(merged$id==merged$id.true),]))/length(rownames(merged[which(!is.na(merged$id)),]))
recall_all <- length(rownames(merged[which(merged$id==merged$id.true),]))/length(rownames(merged))

prec_S <- length(rownames(merged[which(merged$id==merged$id.true & merged$type=="S"),]))/length(rownames(merged[which(merged$type=="S"),]))
recall_S <- length(rownames(merged[which(merged$id==merged$id.true & merged$type.true=="S"),]))/length(rownames(merged[which(merged$type.true=="S"),]))
prec_D <- length(rownames(merged[which(merged$id==merged$id.true & merged$type=="D"),]))/length(rownames(merged[which(merged$type=="D"),]))
recall_D <- length(rownames(merged[which(merged$id==merged$id.true & merged$type.true=="D"),]))/length(rownames(merged[which(merged$type.true=="D"),]))

df <- data.frame(prec=c(prec_all), recall=c(recall_all), prec_S=c(prec_S), recall_S=c(recall_S), prec_D=c(prec_D), recall_D=c(recall_D), llrmean=c(meanllr), prog=c(prog), donor=c(donor), num=c(num))
write.table(df, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

