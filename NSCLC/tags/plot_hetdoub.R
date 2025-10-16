#! /usr/bin/env Rscript
library(ggplot2)
library(ggsci)

args <- commandArgs(trailingOnly=TRUE)

all <- read.table('all_hetdoub_ambsamp.stats')

all$V1 <- gsub("ambsamp_", "", all$V1)
all$V1 <- gsub("_", ".", all$V1, fixed=T)
all$V1 <- as.numeric(all$V1)

all2 <- read.table('all_hetdoub.stats')
all2$V1 <- gsub("s|b|c", "", all2$V1)
all2$V1 <- gsub("_", ".", all2$V1, fixed=T)
#all2$V1 <- as.numeric(all2$V1)
all <- rbind(all, all2)

colnames(all) <- c("num", "prog", "hetdoub")
means <- aggregate(all$hetdoub, by=list(num=all$num, prog=all$prog), FUN=mean)
mins <- aggregate(all$hetdoub, by=list(num=all$num, prog=all$prog), FUN=min)
maxes <- aggregate(all$hetdoub, by=list(num=all$num, prog=all$prog), FUN=max)

colnames(means)[3] <- "mean"
colnames(mins)[3] <- "min"
colnames(maxes)[3] <- "max"

progs_alphabetical <- unique(all[order(all$prog),]$prog)
progs <- progs_alphabetical[grep("CellBouncer", progs_alphabetical, invert=T)]
progs_cbfirst <- c("CellBouncer", progs)
progs_cblast <- c(progs, "CellBouncer")
all$prog <- factor(all$prog, labels=progs_cblast, levels=progs_cblast)
all <- all[order(all$prog),]
#all$prog <- factor(all$prog, labels=progs_cbfirst, levels=progs_cbfirst)
pd <- position_dodge()

pal <- pal_observable()(length(progs_cbfirst))
name2col <- setNames(pal, progs_cbfirst)

plt <- ggplot(all) + 
    geom_violin(aes(x=prog, y=hetdoub, fill=prog, group=prog)) +
    scale_y_continuous("Heterotypic doublet rate", breaks=seq(0,1,0.1)) + 
    theme_bw() + 
    geom_hline(aes(yintercept=0.1572569), lty='dotted') + 
    scale_fill_manual("Program", values=name2col) +
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12, angle = 90, vjust = 0.5, hjust=1), 
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14), 
          axis.title.y=element_text(size=14),
          strip.background=element_blank(),
          strip.text.y=element_text(face='bold', size=14),
          legend.position='top')

ggsave(plt, file="hetdoub_amb.pdf", width=5, height=4)


