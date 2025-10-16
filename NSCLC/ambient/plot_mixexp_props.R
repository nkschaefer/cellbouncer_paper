#! /usr/bin/env Rscript
library(ggplot2)
library(ggthemes)

dat <- read.table('../demux_vcf_sc_1000.contam_prof')
dat$perc <- 0

dat1 <- dat
dat1$experiment <- "donor_1"
dat2 <- dat
dat2$experiment <- "donor_2"
dat3 <- dat
dat3$experiment <- "donor_3"
dat4 <- dat
dat4$experiment <- "donor_4"
dat5 <- dat
dat5$experiment <- "donor_5"
dat6 <- dat
dat6$experiment <- "donor_6"
dat7 <- dat
dat7$experiment <- "donor_7"
dat <- rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7)

datcorr <- dat

bp <- read.table('../demux_vcf_sc_1000.bulkprops')
bp <- bp[,c(1,2)]
bp$perc <- 0
bp1 <- bp
bp1$experiment <- "donor_1"
bp2 <- bp
bp2$experiment <- "donor_2"
bp3 <- bp
bp3$experiment <- "donor_3"
bp4 <- bp
bp4$experiment <- "donor_4"
bp5 <- bp
bp5$experiment <- "donor_5"
bp6 <- bp
bp6$experiment <- "donor_6"
bp7 <- bp
bp7$experiment <- "donor_7"
bp <- rbind(bp1, bp2, bp3, bp4, bp5, bp6, bp7)

cellcounts2props <- function(name){
    dat <- read.table(name)
    s <- dat[which(dat$V3=="S"),]
    d <- dat[which(dat$V3=="D"),]
    d$id1 <- apply(d, 1, function(x){ strsplit(x[2], "+", fixed=T)[[1]][1]})
    d$id2 <- apply(d, 1, function(x){ strsplit(x[2], "+", fixed=T)[[1]][2]})
    s <- s[,c(1,2)]
    d1 <- d[,c(1,5)]
    d2 <- d[,c(1,6)]
    colnames(s) <- c("bc", "id")
    colnames(d1) <- c("bc", "id")
    colnames(d2) <- c("bc", "id")
    s$num <- 2
    d1$num <- 1
    d2$num <- 1
    all <- rbind(s, d1, d2)
    agg <- aggregate(all$num, by=list(id=all$id), FUN=sum)
    agg$x <- agg$x / sum(agg$x)
    colnames(agg)[2] <- "cells"
    return(agg)
}

cells <- cellcounts2props('../demux_vcf_sc_1000.assignments')
cells$perc <- 0
c1 <- cells
c1$experiment <- "donor_1"
c2 <- cells
c2$experiment <- "donor_2"
c3 <- cells
c3$experiment <- "donor_3"
c4 <- cells
c4$experiment <- "donor_4"
c5 <- cells
c5$experiment <- "donor_5"
c6 <- cells
c6$experiment <- "donor_6"
c7 <- cells
c7$experiment <- "donor_7"
cells <- rbind(c1, c2, c3, c4, c5, c6, c7)


ratebase <- read.table('../demux_vcf_sc_1000.contam_rate')
c <- median(ratebase$V2)
rates <- data.frame(rate=c(c,c,c,c,c,c,c), perc=c(0,0,0,0,0,0,0), 
    id=c("donor_1", "donor_2", "donor_3", "donor_4", "donor_5", "donor_6", "donor_7"))

#rates <- NA
#first <- TRUE
first <- FALSE

percs <- c('0.05', '0.1', '0.15', '0.2', '0.25')
for (d in seq(1, 7)){
    for (perc in percs){
        name <- paste('d', d, 'c', perc, '.contam_prof', sep='')
        tab <- read.table(name)
        tab$perc <- as.numeric(perc)
        tab$experiment <- paste('donor_', d, sep='')
        
        tabcorr <- read.table(paste('d', d, 'c', perc, '_corr.contam_prof', sep=''))
        tabcorr$perc <- as.numeric(perc)
        tabcorr$experiment <- paste('donor_', d, sep='')

        tabrate <- read.table(paste('d', d, 'c', perc, '.contam_rate', sep=''))
        row <- data.frame(rate=median(tabrate$V2), perc=as.numeric(perc), id=paste('donor_', d, sep=''))
        
        bpthis <- read.table(paste('d', d, 'c', perc, '.bulkprops', sep=''))
        bpthis <- bpthis[,c(1,2)]
        bpthis$perc <- as.numeric(perc)
        bpthis$experiment <- paste('donor_', d, sep='')
        
        cellsthis <- cellcounts2props(paste('d', d, 'c', perc, '.decontam.assignments', sep=''))
        cellsthis$perc <- as.numeric(perc)
        cellsthis$experiment <- paste('donor_', d, sep='')

        if (first){
            dat <- tab
            rates <- row
            bp <- bpthis
            cells <- cellsthis
            first <- FALSE
        } else{
            dat <- rbind(dat, tab)
            rates <- rbind(rates, row)
            bp <- rbind(bp, bpthis)
            cells <- rbind(cells, cellsthis)
            datcorr <- rbind(datcorr, tabcorr)
        }
    }
}


#dat <- dat[which(dat$V1==dat$id),]
colnames(dat) <- c("id", "p", "p.conc", "perc", "experiment")
dat$expected <- dat$perc/(dat$perc + c)

datcorr <- datcorr[,c(1,2,4,5)]
colnames(datcorr) <- c("id", "p.corr", "perc", "experiment")

colnames(bp) <- c("id", "bulk", "perc", "experiment")
dat <- merge(dat, bp)
dat <- merge(dat, datcorr)
dat <- merge(dat, cells)
write.table(dat, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

library(reshape2)
melted <- melt(dat, id.vars=c('id', 'perc', 'experiment'))
melted <- melted[which(melted$variable != 'p.conc' & melted$variable != "expected"),]
#melted$variable <- factor(melted$variable, labels=c("Bulk", "Cells", "Ambient"), levels=c("bulk", "cells", "p"))
melted$variable <- factor(melted$variable, labels=c("Bulk", "Cells", "Ambient", "Amb_corr"), levels=c("bulk", "cells", "p", "p.corr"))

colnames(melted)[1] <- "Individual"

plt <- ggplot(melted[which(melted$perc < 0.3),]) + 
    geom_bar(aes(x=perc, y=value, fill=Individual), width=0.05, stat='identity') + 
    facet_grid(variable~experiment) + 
    scale_x_continuous("Percent introduced") +
    scale_y_continuous("Fraction") + 
    scale_fill_tableau() +
    theme_bw() + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          strip.background=element_blank(),
          strip.text.x=element_text(size=14, face='bold'),
          strip.text.y=element_text(size=14, face='bold'),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14))
ggsave(plt, file="mixexp_props_true.pdf", width=12, height=6)

colnames(rates)[3] <- "Individual"
plt <- ggplot(rates) +
    geom_abline(aes(slope=1, intercept=0), lty='dotted', lwd=1.5) + 
    geom_point(aes(x=perc, y=rate, colour=Individual), size=2.5) +
    scale_x_continuous("Percent introduced") + 
    scale_y_continuous("Median contam. rate") + 
    scale_colour_tableau() + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12, angle=90, vjust=0.5, hjust=1),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14))
ggsave(plt, file="mixexp_rates.pdf", width=5, height=3)


