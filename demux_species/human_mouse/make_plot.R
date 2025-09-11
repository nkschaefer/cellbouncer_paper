#! /usr/bin/env Rscript
library(ggplot2)
library(ggthemes)

args <- commandArgs(trailingOnly=TRUE)

dat <- read.table('mprofile_run_all.txt')

mins <- aggregate(dat$V2, by=list(prog=dat$V3), FUN='min')
colnames(dat) <- c("mem", "time", "prog")
dat <- merge(dat, mins)
dat$time2 <- dat$time - dat$x

dat2 <- dat

progs <- c("STARsolo", "k52", "k42", "k32", "k22")
dat2$prog <- factor(dat2$prog, labels=progs, levels=progs)
dat2 <- dat2[order(dat2$prog),]

plt <- ggplot(dat2) + 
    geom_point(aes(x=time2/60, y=mem/1024, colour=prog)) + 
    geom_line(aes(x=time2/60, y=mem/1024, colour=prog, group=prog)) + 
    theme_few() + 
    scale_colour_tableau() + 
    scale_x_continuous("Time (minutes)") + 
    scale_y_continuous("Memory (GB)") + 
    theme(axis.text.x=element_text(size=12), 
          axis.title.x=element_text(size=14), 
          axis.text.y=element_text(size=12), 
          axis.title.y=element_text(size=14))

pltname = "mem_run.pdf"
ggsave(plt, file=pltname, width=8, height=5)

