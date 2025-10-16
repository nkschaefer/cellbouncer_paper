#! /usr/bin/env Rscript
library(ggplot2)
library(ggsci)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1){
    write("Please provide stats file", stderr())
    q()
}
all <- read.table(args[1], sep='\t')

all$V1 <- gsub("ambsamp_", "", all$V1)
all$V1 <- gsub("s|b|c", "", all$V1)
all$V1 <- gsub("_", ".", all$V1, fixed=T)
all$V1 <- as.numeric(all$V1)

all$f <- 2*(all$V3*all$V4)/(all$V3+all$V4)

colnames(all)[1] <- "num"

meanp <- aggregate(all$V3, by=list(num=all$num, prog=all$V2), FUN=mean)
minp <- aggregate(all$V3, by=list(num=all$num, prog=all$V2), FUN=min)
maxp <- aggregate(all$V3, by=list(num=all$num, prog=all$V2), FUN=max)

meanr <- aggregate(all$V4, by=list(num=all$num, prog=all$V2), FUN=mean)
minr <- aggregate(all$V4, by=list(num=all$num, prog=all$V2), FUN=min)
maxr <- aggregate(all$V4, by=list(num=all$num, prog=all$V2), FUN=max)

meanf <- aggregate(all$f, by=list(num=all$num, prog=all$V2), FUN=mean)
minf <- aggregate(all$f, by=list(num=all$num, prog=all$V2), FUN=min)
maxf <- aggregate(all$f, by=list(num=all$num, prog=all$V2), FUN=max)

colnames(meanp)[3] <- "mean"
colnames(minp)[3] <- "min"
colnames(maxp)[3] <- "max"
colnames(meanr)[3] <- "mean"
colnames(minr)[3] <- "min"
colnames(maxr)[3] <- "max"
colnames(meanf)[3] <- "mean"
colnames(minf)[3] <- "min"
colnames(maxf)[3] <- "max"

pmerge <- merge(meanp, minp)
pmerge <- merge(pmerge, maxp)
pmerge$stat <- "Precision"

rmerge <- merge(meanr, minr)
rmerge <- merge(rmerge, maxr)
rmerge$stat <- "Recall"

fmerge <- merge(meanf, minf)
fmerge <- merge(fmerge, maxf)
fmerge$stat <- "F"

all <- rbind(pmerge, rmerge, fmerge)

progs_alphabetical <- unique(all[order(all$prog),]$prog)
progs <- progs_alphabetical[grep("CellBouncer", progs_alphabetical, invert=T)]
#progs_cbfirst <- c(c("CellBouncer_filtered", "CellBouncer_unfiltered"), progs)
#progs_cblast <- c(progs, c("CellBouncer_filtered", "CellBouncer_unfiltered"))
progs_cbfirst <- c("CellBouncer", progs)
progs_cblast <- c(progs, "CellBouncer")
all$prog <- factor(all$prog, labels=progs_cblast, levels=progs_cblast)
all <- all[order(all$prog),]
#all$prog <- factor(all$prog, labels=progs_cbfirst, levels=progs_cbfirst)
pd <- position_dodge()

pal <- pal_observable()(length(progs_cbfirst))
name2col <- setNames(pal, progs_cbfirst)

plt <- ggplot(all[which(all$stat=="F"),]) + 
    geom_line(aes(x=num, y=mean, colour=prog, group=prog), lwd=1.5) + 
    geom_segment(aes(x=num, xend=num, y=min, yend=max, colour=prog, group=prog)) + 
    geom_point(aes(x=num, y=mean, colour=prog), size=2) + 
    #geom_errorbar(aes(x=num, ymin=min, ymax=max, colour=prog, group=prog), position=pd) + 
    #geom_pointrange(aes(x=num, y=mean, ymin=min, ymax=max, colour=prog), position=pd) + 
    facet_grid(stat~.) +
    scale_x_continuous("Percent foreground tags removed") + 
    #scale_x_log10("Percent foreground tags removed") + 
    #annotation_logticks(side="b") +
    scale_y_continuous("Mean value", breaks=seq(0,1,0.1)) + 
    theme_bw() + 
    scale_color_manual("Program", values=name2col) +
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),
          strip.background=element_blank(),
          strip.text.y=element_text(face='bold', size=14),
          legend.position='top')

#plt <- ggplot(all) + geom_bar(aes(x=paste(num, prog), fill=prog, y=mean), stat='identity') + 
#    geom_errorbar(aes(x=paste(num, prog), ymin=min, ymax=max, colour=prog)) + 
#    facet_grid(stat~num, scales='free_x') + 
#    theme_bw() + 
#    theme(panel.grid.major=element_blank(), 
#        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
#         axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))

ggsave(plt, file="stats_ambsamp.pdf", width=5, height=3.5)

all <- all[order(all$mean, decreasing=T),]
all <- all[order(all$stat),]
all <- all[order(all$num),]

write.table(all, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

