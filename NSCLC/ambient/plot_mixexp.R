#! /usr/bin/env Rscript
library(ggplot2)

cols <- c("#EF476F", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")
progs_ord <- c("demuxalot", "demuxlet", "vireo", "CellBouncer", "CellBouncer+contam")
names(cols) <- progs_ord

dat <- read.table('mixexp_compiled.txt')
dat <- dat[,c(1,2,3,4,5)]
dat[which(dat$V1=="CellBouncer_contam"),]$V1 <- "CellBouncer+contam"
colnames(dat) <- c("prog", "perc", "donor", "prec", "recall")
dat$f <- 2*(dat$prec*dat$recall)/(dat$prec+dat$recall)

fmean <- aggregate(dat$f, by=list(prog=dat$prog, perc=dat$perc), FUN=mean)
fmin <- aggregate(dat$f, by=list(prog=dat$prog, perc=dat$perc), FUN=min)
fmax <- aggregate(dat$f, by=list(prog=dat$prog, perc=dat$perc), FUN=max)

colnames(fmean)[3] <- "f"
colnames(fmin)[3] <- "fmin"
colnames(fmax)[3] <- "fmax"

merged <- merge(fmean, fmin)
merged <- merge(merged, fmax)

dat2 <- read.table('../sens_spec_nvar.txt')
dat2 <- dat2[which(dat2$V10==500),]
dat2$f <- 2*(dat2$V1*dat2$V2)/(dat2$V1+dat2$V2)
dat2[which(dat2$V8=="CellBouncer_contam"),]$V8 <- "CellBouncer+contam"

addrow <- data.frame(prog=dat2$V8, perc=0, f=dat2$f, fmin=dat2$f, fmax=dat2$f)
merged <- rbind(merged, addrow)

merged$prog <- factor(merged$prog, labels=progs_ord, levels=progs_ord)
merged <- merged[order(merged$prog),]


plt <- ggplot(merged) +
    geom_line(aes(x=perc, y=f, colour=prog, group=prog), lwd=1.5) +
    geom_point(aes(x=perc, y=f, colour=prog, group=prog), size=2.5) + 
    geom_errorbar(aes(x=perc, ymin=fmin, ymax=fmax, colour=prog), width=0.005) +
    scale_colour_manual("Program", values=cols) +
    theme_bw() +
    scale_x_continuous("Percent introduced contamination") +
    scale_y_continuous("F score") +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          legend.title=element_text(size=14),
          legend.text=element_text(size=12),
          legend.position="top")

ggsave(plt, file="mixexp_results.pdf", width=5, height=4)

