#! /usr/bin/env Rscript

# Load data 
dat <- read.table('inter_species_dat.tsv', header=T)
dat2 <- read.table('intra_species_dat.tsv', header=T)

spec <- read.table('../bc2species.tsv', sep='\t')
colnames(spec) <- c("bc", "species")

dat <- dat[,which(colnames(dat) != "species")]
dat2 <- dat2[,which(colnames(dat2) != "species")]

dat <- merge(dat, spec)
dat2 <- merge(dat2, spec)

dat <- dat[which((dat$s1=="Human" & dat$s2=="Chimp" & dat$species=="Human_Chimp") | 
                 (dat$s1=="Chimp" & dat$s2=="Bonobo" & dat$species=="Chimp_Bonobo")),]
dat2 <- dat2[which(dat2$species %in% c("Human_Human", "Chimp_Chimp")),]

c1 <- dat2[which(dat2$id2 < dat2$id1),]$count1
dat2[which(dat2$id2 < dat2$id1),]$count1 <- dat2[which(dat2$id2 < dat2$id1),]$count2
dat2[which(dat2$id2 < dat2$id1),]$count2 <- c1


s11 <- aggregate(dat$count1, by=list(bc=dat$bc, species=dat$species, lib=dat$lib), FUN=sum)
s12 <- aggregate(dat$count2, by=list(bc=dat$bc, species=dat$species, lib=dat$lib), FUN=sum)

s21 <- aggregate(dat2$count1, by=list(bc=dat2$bc, species=dat2$species, lib=dat2$lib), FUN=sum)
s22 <- aggregate(dat2$count2, by=list(bc=dat2$bc, species=dat2$species, lib=dat2$lib), FUN=sum)

colnames(s11)[4] <- "count1"
colnames(s12)[4] <- "count2"
colnames(s21)[4] <- "count1"
colnames(s22)[4] <- "count2"

merge1 <- merge(s11, s12)
merge2 <- merge(s21, s22)

all <- rbind(merge1, merge2)
all$frac <- all$count1/(all$count1+all$count2)

a_all <- read.table('../bc2id.tsv', sep='\t')
colnames(a_all) <- c("bc", "id")
all <- merge(all, a_all)

all$id <- gsub("+", "x", all$id, fixed=T)

n1 <- apply(all[which(all$species=="Human_Chimp"),], 1, function(x){
                strsplit(x[7], 'x')[[1]][1] })
n2 <- apply(all[which(all$species=="Human_Chimp"),], 1, function(x){ 
                strsplit(x[7], 'x')[[1]][2] })
all[which(all$species=="Human_Chimp"),]$id <- paste(n2, n1, sep='x')

counts <- as.data.frame(table(all$id))
colnames(counts) <- c("id", "ncells")
all <- merge(all, counts)

key <- unique(all[,which(colnames(all) %in% c("id", "species", "ncells"))])
key <- key[order(key$ncells),]
key <- key[order(key$species),]
all$id <- factor(all$id, labels=key$id, levels=key$id)

all <- all[order(all$frac),]
all <- all[order(all$id),]
all$num <- seq(1, length(rownames(all)))
minnum <- aggregate(all$num, by=list(id=all$id), FUN=min)
colnames(minnum)[2] <- "min"
all <- merge(all, minnum)
all$id1 <- apply(all, 1, function(x){ strsplit(x[1], 'x')[[1]][1] })
all$id2 <- apply(all, 1, function(x){ strsplit(x[1], 'x')[[1]][2] })

chimps <- c("C3624", "C3651", "C40670", "C8861", "C6007B", "C40280")
humans <- c("H20961", "H21792" ,"H23555", "H28126", "H29089")
bonobo <- c("CongoA4B")
all$species1 <- ""
all$species2 <- ""
all[which(all$id1 %in% chimps),]$species1 <- "Chimp"
all[which(all$id1 %in% humans),]$species1 <- "Human"
#all[which(all$id1 %in% bonobo),]$species1 <- "Bonobo"
all[which(all$id2 %in% chimps),]$species2 <- "Chimp"
all[which(all$id2 %in% humans),]$species2 <- "Human"
all[which(all$id2 %in% bonobo),]$species2 <- "Bonobo"
annodat <- unique(all[,c(1,3,11,12,13,14)])

library(ggplot2)
library(viridis)
library(pals)

cols <- c(Chimp="#EF9943", Bonobo="#7C3626", Human="#7EBDC2")

plt <- ggplot(all) +
    geom_tile(aes(x=num-min, y=id, fill=frac)) + 
    geom_text(data=annodat, aes(x=-950, y=id, label=id1, colour=species1), 
        show.legend=FALSE, hjust='left', fontface='bold', size=11, size.unit='pt', family='Helvetica') +
    geom_text(data=annodat, aes(x=-550, y=id, label=id2, colour=species2), 
        show.legend=FALSE, hjust='left', fontface='bold', size=11, size.unit='pt', family='Helvetica') + 
    theme_bw() + 
    scale_x_continuous("Cells") + 
    scale_y_discrete("Contributor cell lines") + 
    #scale_fill_viridis("Percent MT reads from cell line 1") + 
    scale_fill_gradientn("Percent MT reads from line 1", colours=ocean.delta(100), limits=c(0,1)) + 
    #scale_fill_distiller(palette="ocean.delta") +
    #scale_fill_viridis("Percent MT reads from line 1", option='mako') + 
    #scale_fill_gradientn("Percent MT reads from line 1", colours=ocean.curl(1000), limits=c(0,1)) +
    scale_colour_manual("Species", values=cols) +  
    facet_grid(species~., scales='free_y', space='free_y') + 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_blank(),
          #axis.text.y=element_text(size=12),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          strip.background=element_blank(),
          strip.text=element_blank(),
          panel.border=element_blank(),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14))

ggsave(plt, file="mito_fracs.pdf", width=8, height=5, bg='white')

