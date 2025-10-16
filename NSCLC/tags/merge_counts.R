#! /usr/bin/env Rscript
library(reshape2)

combined <- NULL
for (i in seq(1, 12)){
    countsn <- paste("ALL", i, 'counts', sep='.')
    tab <- read.table(countsn, sep="\t", header=T)
    if (i == 1){
        combined <- tab
    }
    else{
        combined <- rbind(combined, tab)
    }
}

melted <- melt(combined, id.vars=c("barcode"))

agg <- aggregate(melted$value, by=list(barcode=melted$barcode, variable=melted$variable), FUN=sum)

casted <- dcast(agg, barcode ~ variable, value.var='x')
write.table(casted, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)


