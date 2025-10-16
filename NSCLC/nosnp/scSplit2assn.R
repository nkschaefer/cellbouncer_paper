#! /usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("Please provide output dir", stderr())
    q()
}

lastchr <- substr(args[1], nchar(args[1]), nchar(args[1]))
if (lastchr == '/'){
    args[1] <- substr(args[1], 1, nchar(args[1]) - 1)
}

fn1 <- paste(args[1], '/', 'scSplit_result.csv', sep='')
fn2 <- paste(args[1], '/', 'scSplit_P_s_c.csv', sep='')

tab1 <- read.table(fn1, sep='\t', header=T)
tab2 <- read.table(fn2, sep=',', header=T)

tab1$type <- "S"
tab1[grep("SNG", tab1$Cluster, invert=T),]$type <- "D"
tab1$Cluster <- gsub("SNG-", "", gsub("DBL-", "", tab1$Cluster))

tab2$top <- apply(tab2, 1, function(x){ colnames(tab2)[ which.max(as.numeric(x[2:length(x)])) + 1 ]})
tab2$second <- apply(tab2, 1, function(x){ colnames(tab2)[ order(as.numeric(x[2:(length(x)-1)]), decreasing=T)[2] + 1 ]})
tab2$third <- apply(tab2, 1, function(x){ colnames(tab2)[ order(as.numeric(x[2:(length(x)-2)]), decreasing=T)[3] + 1 ]})

tab2 <- tab2[,which(colnames(tab2) %in% c("X", "top", "second", "third"))]
colnames(tab2) <- c("Barcode", "top", "second", "third")
tab2$top <- gsub("X", "", tab2$top)
tab2$second <- gsub("X", "", tab2$second)
tab2$third <- gsub("X", "", tab2$third)
tab2 <- tab2[which(tab2$Barcode %in% tab1[which(tab1$type=="D"),]$Barcode),]

# Top choice for doublets is always a distinct cluster for doublets - 
# assume next 2 top choices represent singlet components
tab2$id1 <- tab2$second
tab2$id2 <- tab2$third
tab2[which(tab2$third < tab2$second),]$id1 <- tab2[which(tab2$third < tab2$second),]$third
tab2[which(tab2$third < tab2$second),]$id2 <- tab2[which(tab2$third < tab2$second),]$second

tab2$id <- paste(tab2$id1, tab2$id2, sep='+')
tab2 <- tab2[,which(colnames(tab2) %in% c("Barcode", "id"))]
colnames(tab2)[2] <- "Cluster"
tab2$type <- "D"

tab1 <- tab1[which(tab1$type=="S"),]

tab <- rbind(tab1, tab2)
tab$ll <- 100

write.table(tab, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


