#! /usr/bin/env Rscript
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1){
    write("ERROR: provide GnomAD tsv", stderr())
    q()
}
gnomad <- read.table(args[1], sep='\t', header=T)
library(ggplot2)

mut <- "M-10478-C-T"
mut2 <- "M-10478-C-A"

gnomad <- gnomad[which(gnomad$Homoplasmic.Allele.Frequency > 0),]

high <- c("start lost", "stop lost", "frameshift")
moderate <- c("inframe deletion", "missense")
low <- c("incomplete terminal codon", "start retained", "stop retained", "synonymous", "coding sequence")
modifier <- c("non coding transcript exon")

impact <- data.frame(VEP.Annotation=c(high, moderate, low, modifier), 
                     VEP_IMPACT=c(rep("HIGH", length(high)),
                                  rep("MODERATE", length(moderate)),
                                  rep("LOW", length(low)),
                                  rep("MODIFIER", length(modifier))))

gnomad <- merge(gnomad, impact, all.x=TRUE)

gnomad[which(gnomad$VEP_IMPACT == "MODERATE"),]$VEP_IMPACT <- "MODERATE_HIGH"
gnomad[which(gnomad$VEP_IMPACT == "HIGH"),]$VEP_IMPACT <- "MODERATE_HIGH"
vep <- gnomad[which(gnomad$VEP_IMPACT %in% c("LOW", "MODERATE_HIGH")),]

gnomad$Clinvar.type <- ""
gnomad[which(gnomad$ClinVar.Germline.Classification %in% c("association not found", "Benign/Likely benign", "Benign", "Likely benign")),]$Clinvar.type <- "Benign"
gnomad[which(gnomad$ClinVar.Germline.Classification %in% c("Affects", "Likely pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic", "Likely pathogenic; drug response", "Pathogenic; drug response")),]$Clinvar.type <- "Pathogenic"

clinvar <- gnomad[which(gnomad$Clinvar.type != ""),]

vep <- vep[,which(colnames(vep) %in% c("Homoplasmic.Allele.Frequency", "VEP_IMPACT"))]
clinvar <- clinvar[,which(colnames(clinvar) %in% c("Homoplasmic.Allele.Frequency", "Clinvar.type"))]
colnames(vep) <- c("MAF", "Class")
colnames(clinvar) <- c("MAF", "Class")
vep$type <- "VEP"
clinvar$type <- "Clinvar"

cols <- c(Benign="#B37BA4", Pathogenic="#FF9B71", LOW="#B37BA4", MODERATE_HIGH="#FF9B71")

plt <- ggplot(rbind(vep, clinvar)) + 
    geom_violin(aes(x=Class, y=MAF, fill=Class), show.legend=FALSE, lwd=1.05, draw_quantiles=c(0.5)) + 
    geom_hline(aes(yintercept=8.860064e-05), lty='dotted', lwd=1.1) + 
    facet_grid(.~type, scales='free_x') +
    scale_y_log10("Minor allele frequency") + 
    scale_x_discrete("") + 
    annotation_logticks(side="l") +  
    scale_fill_manual(values=cols) + 
    theme_bw() +
    ggtitle("chrM:10478:C>T") + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          strip.background=element_blank(),
          plot.title=element_text(size=14, face='bold'),
          strip.text.x=element_text(size=12, face='bold'))

ggsave(plt, file="gnomadplot.pdf", width=3.5, height=4.5) 

