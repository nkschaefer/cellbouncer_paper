#! /usr/bin/env Rscript

dat <- read.table('mito_expr_dat.tsv', header=T)
strand <- read.table('mt_strand.txt')
colnames(strand) <- c("gene", "strand")

dat <- dat[grep("^MT-", dat$gene),]
dat <- merge(dat, strand)
dat <- dat[which((dat$id1 == "Human_HC" & dat$id2 == "Chimp_HC") | (dat$id1=="Human" & dat$id2 == "Chimp") | 
                 (dat$id1=="Chimp_HC" & dat$id2 == "Human_HC")),]
library(reshape2)
casted <- dcast(dat, gene + strand ~ id1 + id2, value.var='lfc')
casted$MT <- casted$Human_HC_Chimp_HC
casted$nuclear <- casted$Human_Chimp + casted$Chimp_HC_Human_HC

print("MT KS test")
print(ks.test(casted[which(casted$strand=="Heavy"),]$MT, casted[which(casted$strand=="Light"),]$MT))
print("MT wilcox test")
print(wilcox.test(casted[which(casted$strand=="Heavy"),]$MT, casted[which(casted$strand=="Light"),]$MT))
print("nuclear KS test")
print(ks.test(casted[which(casted$strand=="Heavy"),]$nuclear, casted[which(casted$strand=="Light"),]$nuclear))
print("nuclear wilcox test")
print(wilcox.test(casted[which(casted$strand=="Heavy"),]$nuclear, casted[which(casted$strand=="Light"),]$nuclear))


melted <- melt(casted[,which(colnames(casted) %in% c("gene", "strand", "MT", "nuclear"))], id.vars=c("gene", "strand"))
library(ggplot2)
cols <- c(Heavy="#B37BA4", Light="#FF9B71")
plt <- ggplot(melted) + 
    geom_violin(aes(x=strand, y=value, fill=strand), lwd=1.05, show.legend=FALSE, draw_quantiles=c(0.5)) + 
    facet_grid(.~variable) + 
    theme_bw() + 
    scale_x_discrete("Strand") + 
    scale_y_continuous("Log2 Fold Change") + 
    scale_fill_manual(values=cols) + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          strip.text.y=element_text(size=14, face='bold'),
          strip.text.x=element_text(size=14, face='bold'),
          strip.background=element_blank())
ggsave(plt, file="heavy_light.pdf", width=3.5, height=3)

