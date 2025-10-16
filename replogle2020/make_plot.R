#! /usr/bin/env Rscript

library(ggplot2)
library(viridis)
library(reshape2)

cbfull <- read.table('replogle_combined.table', header=T)
cbfull$count <- rowSums(cbfull[,-c(1)])
cbfull <- cbfull[,which(colnames(cbfull) %in% c("barcode", "count"))]

cb50 <- read.table('replogle_combined_ds0.5.table', header=T)
cb50$count <- rowSums(cb50[,-c(1)])
cb50 <- cb50[,which(colnames(cb50) %in% c("barcode", "count"))]

cb25 <- read.table('replogle_combined_ds0.25.table', header=T)
cb25$count <- rowSums(cb25[,-c(1)])
cb25 <- cb25[,which(colnames(cb25) %in% c("barcode", "count"))]

cb10 <- read.table('replogle_combined_ds0.1.table', header=T)
cb10$count <- rowSums(cb10[,-c(1)])
cb10 <- cb10[,which(colnames(cb10) %in% c("barcode", "count"))]

cb5 <- read.table('replogle_combined_ds0.05.table', header=T)
cb5$count <- rowSums(cb5[,-c(1)])
cb5 <- cb5[,which(colnames(cb5) %in% c("barcode", "count"))]

cb1 <- read.table('replogle_combined_ds0.01.table', header=T)
cb1$count <- rowSums(cb1[,-c(1)])
cb1 <- cb1[,which(colnames(cb1) %in% c("barcode", "count"))]

pgfull <- read.table('replogle_combined.pg.table', header=T)
pgfull$count <- rowSums(pgfull[,-c(1)])
pgfull <- pgfull[,which(colnames(pgfull) %in% c("barcode", "count"))]

pg50 <- read.table('replogle_combined_ds0.5.pg.table', header=T)
pg50$count <- rowSums(pg50[,-c(1)])
pg50 <- pg50[,which(colnames(pg50) %in% c("barcode", "count"))]

pg25 <- read.table('replogle_combined_ds0.25.pg.table', header=T)
pg25$count <- rowSums(pg25[,-c(1)])
pg25 <- pg25[,which(colnames(pg25) %in% c("barcode", "count"))]

pg10 <- read.table('replogle_combined_ds0.1.pg.table', header=T)
pg10$count <- rowSums(pg10[,-c(1)])
pg10 <- pg10[,which(colnames(pg10) %in% c("barcode", "count"))]

pg5 <- read.table('replogle_combined_ds0.05.pg.table', header=T)
pg5$count <- rowSums(pg5[,-c(1)])
pg5 <- pg5[,which(colnames(pg5) %in% c("barcode", "count"))]

pg1 <- read.table('replogle_combined_ds0.01.pg.table', header=T)
pg1$count <- rowSums(pg1[,-c(1)])
pg1 <- pg1[,which(colnames(pg1) %in% c("barcode", "count"))]


colnames(cbfull)[2] <- "CellBouncer"
colnames(cb50)[2] <- "CellBouncer"
colnames(cb25)[2] <- "CellBouncer"
colnames(cb10)[2] <- "CellBouncer"
colnames(cb5)[2] <- "CellBouncer"
colnames(cb1)[2] <- "CellBouncer"

colnames(pgfull)[2] <- "Mixture_model"
colnames(pg50)[2] <- "Mixture_model"
colnames(pg25)[2] <- "Mixture_model"
colnames(pg10)[2] <- "Mixture_model"
colnames(pg5)[2] <- "Mixture_model"
colnames(pg1)[2] <- "Mixture_model"

merge2 <- function(cb, pg){
    merged <- merge(cb, pg, all.x=TRUE, all.y=TRUE)
    if (length(rownames(merged[which(is.na(merged$CellBouncer)),])) > 0){
        merged[which(is.na(merged$CellBouncer)),]$CellBouncer <- 0
    }
    if (length(rownames(merged[which(is.na(merged$Mixture_model)),])) > 0){
        merged[which(is.na(merged$Mixture_model)),]$Mixture_model <- 0
    }
    return(merged)
}

countfull <- function(){
    c <- 0
    nr <- 0
    for (i in seq(1,6)){
        name <- paste('replogle_split_', i, '.counts', sep='')
        counts <- read.table(name, header=T)
        rs <- rowSums(counts[,-c(1)])
        c <- c + sum(rs)
        nr <- nr + length(rs)
    }
    return(c/nr)
}

full <- merge2(cbfull, pgfull)

ds50 <- merge2(cb50, pg50)
ds25 <- merge2(cb25, pg25)
ds10 <- merge2(cb10, pg10)
ds5 <- merge2(cb5, pg5)
ds1 <- merge2(cb1, pg1)

cfull <- countfull()

cfull <- 999

full$mean <- cfull
ds50$mean <- cfull*0.5
ds25$mean <- cfull*0.25
ds10$mean <- cfull*0.1
ds5$mean <- cfull*0.05
ds1$mean <- cfull*0.01

all <- rbind(full, ds50, ds25, ds10, ds5, ds1)

all$mean <- factor(round(all$mean))

# Count of 6 is 99.9th percentile in CellBouncer
alltrunc <- all[which(all$CellBouncer < 6 & all$Mixture_model < 6),]
alltrunc$count <- 1
agg <- aggregate(alltrunc$count, by=list(CellBouncer=alltrunc$CellBouncer, Mixture_model=alltrunc$Mixture_model, mean=alltrunc$mean), FUN=sum)
colnames(agg)[4] <- "Cells"
colnames(agg)[3] <- "Reads_per_cell"

agg$Cells <- agg$Cells / 1000

for (i in seq(0,5)){
    for (j in seq(0,5)){
        for (rpc in unique(agg$Reads_per_cell)){
            if (length(rownames(agg[which(agg$CellBouncer == i & agg$Mixture_model == j & agg$Reads_per_cell == rpc),])) == 0){
                row <- data.frame(CellBouncer=i, Mixture_model=j, Reads_per_cell=rpc, Cells=0)
                agg <- rbind(agg, row)
            }
        }
    }
}
plt <- ggplot(agg) + geom_tile(aes(x=Mixture_model, y=CellBouncer, fill=Cells)) + 
    theme_bw() + 
    facet_wrap(.~Reads_per_cell, ncol=2) + 
    scale_fill_viridis(option="mako") + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.title.y=element_text(size=14), 
          axis.title.x=element_text(size=14),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.position='top',
          strip.background=element_blank(),
          strip.text.x=element_text(size=14, face='bold'))

ggsave(plt, file="replogle_count_plot.pdf", width=3, height=5.5)

