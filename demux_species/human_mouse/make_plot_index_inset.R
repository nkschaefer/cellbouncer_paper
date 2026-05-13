#! /usr/bin/env Rscript
library(ggplot2)
library(ggthemes)

dat <- read.table('mprofile_index_all.txt')
mins <- aggregate(dat$V2, by=list(prog=dat$V3), FUN='min')
colnames(dat) <- c("mem", "time", "prog")
dat <- merge(dat, mins)
dat$time2 <- dat$time - dat$x

progs <- c("STARsolo", "k52", "k42", "k32", "k22")
dat$prog <- factor(dat$prog, labels=progs, levels=progs)
dat <- dat[order(dat$prog),]

dat <- dat[which(dat$prog != "STARsolo"),]

plt <- ggplot(dat) + 
    geom_point(aes(x=time2/60, y=mem/1024, colour=prog), show.legend=FALSE) + 
    geom_line(aes(x=time2/60, y=mem/1024, colour=prog, group=prog), show.legend=FALSE) + 
    theme_few() + 
    scale_colour_tableau() + 
    scale_x_continuous("Time (minutes)") + 
    scale_y_continuous("Memory (GB)") + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.x=element_blank(), 
          axis.text.y=element_text(size=12), 
          axis.title.y=element_text(size=14))

ggsave(plt, file="mem_index_inset.pdf", width=1.5, height=3)

