#! /usr/bin/env Rscript
library(ggplot2)
library(reshape2)

dat <- read.table('flank_corr_dat.tsv', header=T)

dat <- melt(dat, id.vars=c("gene1", "gene2"))
colnames(dat)[3] <- "type"
dat <- dat[which(! dat$type %in% c("end2", "start", "end")),]

dat$MT <- ""
dat[which(dat$type %in% c("Human", "Human_HC")),]$MT <- "Human"
dat[which(dat$type %in% c("Chimp", "Chimp_HC", "Chimp_CB")),]$MT <- "Chimp"
dat[which(dat$type == "Human_Chimp"),]$MT <- "Human_Chimp"
dat[which(dat$type == "Chimp_Bonobo"),]$MT <- "Chimp_Bonobo"
dat$Genome <- ""
dat[which(dat$type=="Human"),]$Genome <- "Human"
dat[which(dat$type=="Chimp"),]$Genome <- "Chimp"
dat[which(dat$type %in% c("Human_HC", "Chimp_HC", "Human_Chimp")),]$Genome <- "Human_Chimp"
dat[which(dat$type %in% c("Chimp_Bonobo", "Chimp_CB", "Bonob_CB")),]$Genome <- "Chimp_Bonobo"

anno_dat <- unique(dat[,which(colnames(dat) %in% c("type", "MT", "Genome"))])

ord1 <- c("Chimp_Bonobo", "Bonobo_CB", "Chimp_CB", "Human_Chimp", "Human_HC", "Chimp_HC", "Human", "Chimp")
ord2 <- c("Chimp_Bonobo", "Bonobo_CB", "Chimp_CB", "Human_Chimp", "Human_HC", "Chimp_HC", "Human", "Chimp")
dat$type <- factor(dat$type, labels=ord2, levels=ord1)
anno_dat$type <- factor(anno_dat$type, labels=ord2, levels=ord1)
anno_dat$MT <- factor(anno_dat$MT, labels=ord2, levels=ord1)
anno_dat$Genome <- factor(anno_dat$Genome, labels=ord2, levels=ord1)
cols1=c(Human="#7EBDC2", Chimp="#EF9943", Human_Chimp="#F3DFA2", Chimp_Bonobo="#BB4430")
colsA <- c("Human", "Chimp", "Human_Chimp", "Chimp_Bonobo")
colsB <- c("#7EBDC2", "#EF9943", "#2b4162", "#a63a50")
colsC <- setNames(colsB, colsA)

library(reshape2)
dat <- read.table('corrmap_new3.txt', sep='\t')
colnames(dat) <- c("gene1", "gene2", "pos", "value", "weighted", "min", "type")

add_names <- function(dat){
    dat$MT <- ""
    dat[which(dat$type %in% c("Human", "Human_HC")),]$MT <- "Human"
    dat[which(dat$type %in% c("Chimp", "Chimp_HC")),]$MT <- "Chimp"
    dat[which(dat$type == "Human_Chimp"),]$MT <- "Human_Chimp"
    dat[which(dat$type == "Chimp_Bonobo"),]$MT <- "Chimp_Bonobo"
    dat$Genome <- ""
    dat[which(dat$type=="Human"),]$Genome <- "Human"
    dat[which(dat$type=="Chimp"),]$Genome <- "Chimp"
    dat[which(dat$type %in% c("Human_HC", "Chimp_HC", "Human_Chimp")),]$Genome <- "Human_Chimp"
    dat[which(dat$type=="Chimp_Bonobo"),]$Genome <- "Chimp_Bonobo"
    return(dat)
}

dat$type <- factor(dat$type, labels=ord2, levels=ord1)
casted <- dcast(dat, gene1 + gene2 ~ type, value.var='value')

melt1 <- melt(casted, id.vars=c("gene1", "gene2", "Chimp"))
melt2 <- melt(casted, id.vars=c("gene1", "gene2", "Human"))

melt1 <- melt1[which(melt1$variable %in% c("Chimp_HC", "Chimp_Bonobo", "Human_Chimp")),]
melt2 <- melt2[which(melt2$variable %in% c("Human_HC", "Human_Chimp")),]

melt1$variable <- as.character(melt1$variable)
melt2$variable <- as.character(melt2$variable)

melt1$MT <- melt1$variable
melt1[which(melt1$variable=="Chimp_HC"),]$MT <- "Chimp"
melt1[which(melt1$variable=="Chimp_HC"),]$variable <- "Human_Chimp"
melt2$MT <- melt2$variable
melt2[which(melt2$variable=="Human_HC"),]$MT <- "Human"
melt2[which(melt2$variable=="Human_HC"),]$variable <- "Human_Chimp"

plt <- ggplot(melt1) + 
    geom_abline(aes(slope=1, intercept=0), colour="lightgrey") + 
    geom_point(aes(x=Chimp, y=value, colour=MT), show.legend=FALSE, size=5) + 
    geom_point(data=melt1[which(melt1$variable != melt1$MT),], aes(x=Chimp, y=value, colour=variable), show.legend=FALSE, size=2.5) + 
    theme_bw() + 
    scale_colour_manual(values=colsC) + 
    scale_fill_manual(values=colsC) + 
    scale_x_continuous("Gene expr correlation, Chimp") + 
    scale_y_continuous("Gene expr correlation, Other") + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.line.x=element_line(),
          axis.line.y=element_line(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14)) 
ggsave(plt, file="MT_corr_chimp.pdf", width=4, height=4)

plt <- ggplot(melt2) + 
    geom_abline(aes(slope=1, intercept=0), colour="lightgrey") + 
    geom_point(aes(x=Human, y=value, colour=MT), show.legend=FALSE, size=5) + 
    geom_point(data=melt2[which(melt2$variable != melt2$MT),], aes(x=Human, y=value, colour=variable), show.legend=FALSE, size=2.5) + 
    theme_bw() + 
    scale_colour_manual(values=colsC) + 
    scale_fill_manual(values=colsC) + 
    scale_x_continuous("Gene expr correlation, Human") + 
    scale_y_continuous("Gene expr correlation, Other") + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.line.x=element_line(),
          axis.line.y=element_line(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14)) 
ggsave(plt, file="MT_corr_human.pdf", width=4, height=4)

