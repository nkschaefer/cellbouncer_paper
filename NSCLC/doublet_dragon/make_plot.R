#! /usr/bin/env Rscript

vcf <- read.table('demux_vcf_dd.dd.all')
mt <- read.table('demux_mt_dd.dd.all')
tags <- read.table('demux_tags_dd.dd.all')

bp <- read.table('demux_vcf_donor_all_500k.bulkprops')
bp <- bp[,c(1,2)]
colnames(bp) <- c("id", "value")
bp$type <- "bulkprops"

vcfp <- vcf[which(vcf$V1=="group_0"),c(2,3)]
mtp <- mt[which(mt$V1=="group_0"),c(2,3)]
tagsp <- tags[which(tags$V1=="group_0"),c(2,3)]

colnames(vcfp) <- c("id", "value")
colnames(mtp) <- c("id", "value")
colnames(tagsp) <- c("id", "value")

vcfp$type <- "demux_vcf"
mtp$type <- "demux_mt"
tagsp$type <- "demux_tags"

subs_ids <- function(x){
    return(gsub("CMO308", "donor_7",
        gsub("CMO307", "donor_6",
            gsub("CMO306", "donor_5",
                gsub("CMO304", "donor_4",
                    gsub("CMO303", "donor_3",
                        gsub("CMO302", "donor_2",
                            gsub("CMO301", "donor_1", x))))))))
}

mtp$id <- subs_ids(mtp$id)
tagsp$id <- subs_ids(tagsp$id)
props <- rbind(bp, vcfp, mtp, tagsp)

props$data <- "Proportions"

vcfd <- data.frame(id="doublet_rate", value=vcf[which(vcf$V2=="doublet_rate"),]$V3, type="demux_vcf")
mtd <- data.frame(id="doublet_rate", value=mt[which(mt$V2=="doublet_rate"),]$V3, type="demux_mt")
tagsd <- data.frame(id="doublet_rate", value=tags[which(tags$V2=="doublet_rate"),]$V3, type="demux_tags")

rates <- rbind(vcfd, mtd, tagsd)
rates$data <- "Doublet_rate"

vcfl <- data.frame(id="log_likelihood", value=vcf[which(vcf$V1=="data_set"),]$V3, type="demux_vcf")
mtl <- data.frame(id="log_likelihood", value=mt[which(mt$V1=="data_set"),]$V3, type="demux_mt")
tagsl <- data.frame(id="log_likelihood", value=tags[which(tags$V1=="data_set"),]$V3, type="demux_tags")

tags <- rbind(vcfl, mtl, tagsl)
tags$data <- "Log_likelihood"

all <- rbind(props, rates, tags)
drline <- data.frame(value=0.187848, data="Doublet_rate")

library(ggplot2)
library(ggthemes)

progs <- c("demux_mt", "demux_tags", "demux_vcf", "bulkprops")
all$type <- factor(all$type, labels=progs, levels=progs)

dats <- c("Proportions", "Doublet_rate", "Log_likelihood")
all$data2 <- factor(all$data, labels=dats, levels=dats)
drline$data2 <- factor(drline$data, labels=dats, levels=dats)

plt <- ggplot(all) + 
    geom_bar(aes(y=type, x=value, fill=id), stat='identity', show.legend=FALSE) + 
    geom_vline(data=drline, aes(xintercept=value), lwd=1.5, lty='dotted') +     
    facet_grid(.~data2, scales='free_x') + 
    theme_bw() + 
    scale_fill_tableau() + 
    theme(axis.text.x=element_text(size=12, angle = 90, vjust = 0.5, hjust=1),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=12),
          axis.title.y=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          strip.background=element_blank(),
          strip.text.x=element_text(size=12, face='bold'))

ggsave(plt, file="doub_dragon_props.pdf", width=5, height=2)



