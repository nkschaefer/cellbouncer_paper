#! /usr/bin/env Rscript

f <- file('stdin')
dat <- read.table(f)
dat <- dat[,c(4,5,6,7,8,12)]
colnames(dat) <- c("bc", "s1", "s2", "count1", "count2", "gene")
agg1 <- aggregate(dat$count1, by=list(bc=dat$bc, s1=dat$s1, s2=dat$s2, gene=dat$gene), FUN=sum)
agg2 <- aggregate(dat$count2, by=list(bc=dat$bc, s1=dat$s1, s2=dat$s2, gene=dat$gene), FUN=sum)
colnames(agg1)[5] <- "count1"
colnames(agg2)[5] <- "count2"

merged <- merge(agg1, agg2)
write.table(merged, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

