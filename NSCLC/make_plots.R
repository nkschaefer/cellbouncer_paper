#! /usr/bin/env Rscript
library(ggplot2)

cols <- c("#EF476F", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")
progs_ord <- c("demuxalot", "demuxlet", "vireo", "CellBouncer", "CellBouncer+contam")
names(cols) <- progs_ord

dat <- read.table('sens_spec_nvar.txt')
dat[which(dat$V8=="CellBouncer_contam"),]$V8 <- "CellBouncer+contam"

dat <- dat[which(dat$V8 != "CellBouncer+contam"),]
dat$V8 <- factor(dat$V8, labels=progs_ord, levels=progs_ord)
dat <- dat[order(dat$V8),]

dat$f <- 2*(dat$V1*dat$V2)/(dat$V1+dat$V2)

plt <- ggplot(dat) + geom_line(aes(x=V10, y=f, colour=V8, group=V8), lwd=1.5) + 
    geom_point(aes(x=V10, y=f, colour=V8), size=2.5) + 
    scale_colour_manual("Program", values=cols) + 
    theme_bw() +
    scale_x_continuous("SNPs in panel (thousands)") + 
    scale_y_continuous("F score", limits=c(0,1)) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))


ggsave(plt, file="sens_spec_nvar.pdf", width=6, height=4)

tdat <- read.table('time_nvar.txt')
tdat <- tdat[which(tdat$V2 != "CellBouncer_contam"),]
tdat$V2 <- factor(tdat$V2, labels=progs_ord, levels=progs_ord)
tdat <- tdat[order(tdat$V2),]

tplt <- ggplot(tdat) + geom_line(aes(x=V4, y=V1/60/60, colour=V2, group=V2), lwd=1.5) +
    geom_point(aes(x=V4, y=V1/60/60, colour=V2), size=2.5) + 
    scale_colour_manual("Program", values=cols) +
    theme_bw() +
    scale_x_continuous("SNPs in panel (thousands)") +
    scale_y_continuous("Execution time (hours)") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))

ggsave(tplt, file="time_nvar.pdf", width=6, height=4)

dat2 <- read.table('sens_spec_downsample.txt')
dat2$V10 <- as.numeric(gsub("decontam.", "", gsub("ds", "", gsub(".assignments", "", dat2$V9))))
dat2$V10 <- dat2$V10*50
dat2$f <- 2*(dat2$V1*dat2$V2)/(dat2$V1+dat2$V2)
dat2b <- dat[which(dat$V10==500),]
dat2$V8 <- as.character(dat2$V8)
dat2b$V8 <- as.character(dat2b$V8)
dat2b$V10 <- 50

dat2 <- rbind(dat2, dat2b)

dat2[which(dat2$V8=="CellBouncer_contam"),]$V8 <- "CellBouncer+contam"

dat2 <- dat2[which(dat2$V8 != "CellBouncer+contam"),]
dat2$V8 <- factor(dat2$V8, labels=progs_ord, levels=progs_ord)
dat2 <- dat2[order(dat2$V8),]


plt <- ggplot(dat2) + geom_line(aes(x=V10, y=f, colour=V8, group=V8), lwd=1.5) + 
    geom_point(aes(x=V10, y=f, colour=V8), size=2.5) + 
    scale_colour_manual("Program", values=cols) + 
    theme_bw() +
    scale_x_continuous("Expected reads per cell (thousands)") + 
    scale_y_continuous("F score", limits=c(0,1)) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))

ggsave(plt, file="sens_spec_downsample.pdf", width=6, height=4)

tdat2 <- read.table('time_downsample.txt')
tdat2 <- tdat2[which(tdat2$V2 != "CellBouncer_contam"),]
tdat2b <- tdat[which(tdat$V4==500),]
tdat2b$V2 <- as.character(tdat2b$V2)
tdat2b$V4 = 1
tdat2 <- rbind(tdat2, tdat2b)
tdat2$V2 <- factor(tdat2$V2, labels=progs_ord, levels=progs_ord)
tdat2 <- tdat2[order(tdat2$V2),]

tplt <- ggplot(tdat2) + geom_line(aes(x=V4*50, y=V1/60/60, colour=V2, group=V2), lwd=1.5) +
    geom_point(aes(x=V4*50, y=V1/60/60, colour=V2), size=2.5) + 
    scale_colour_manual("Program", values=cols) +
    theme_bw() +
    scale_x_continuous("Expected reads per cell (thousands)") +
    scale_y_continuous("Execution time (hours)") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))

ggsave(tplt, file="time_downsample.pdf", width=6, height=4)


dat3 <- read.table('sens_spec_ncell.txt')
dat3$f <- 2*(dat3$V1*dat3$V2)/(dat3$V1+dat3$V2)
dat3b <- dat[which(dat$V10==500),]
dat3$V8 <- as.character(dat3$V8)
dat3b$V8 <- as.character(dat3b$V8)
dat3b$V10 <- 37000

dat3 <- rbind(dat3, dat3b)

dat3[which(dat3$V8=="CellBouncer_contam"),]$V8 <- "CellBouncer+contam"

dat3 <- dat3[which(dat3$V8 != "CellBouncer+contam"),]
dat3$V8 <- factor(dat3$V8, labels=progs_ord, levels=progs_ord)
dat3 <- dat3[order(dat3$V8),]


plt <- ggplot(dat3) + geom_line(aes(x=V10/1000, y=f, colour=V8, group=V8), lwd=1.5) + 
    geom_point(aes(x=V10/1000, y=f, colour=V8), size=2.5) + 
    scale_colour_manual("Program", values=cols) + 
    theme_bw() +
    scale_x_continuous("Number of cells (thousands)") + 
    scale_y_continuous("F score", limits=c(0,1)) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))

ggsave(plt, file="sens_spec_ncell.pdf", width=6, height=4)

tdat3 <- read.table('time_ncell.txt')
tdat3 <- tdat3[which(tdat3$V2 != "CellBouncer_contam"),]
tdat3b <- tdat[which(tdat$V4==500),]
tdat3b$V2 <- as.character(tdat3b$V2)
tdat3b$V4 = 37000
tdat3 <- rbind(tdat3, tdat3b)
tdat3$V2 <- factor(tdat3$V2, labels=progs_ord, levels=progs_ord)
tdat3 <- tdat3[order(tdat3$V2),]

tplt <- ggplot(tdat3) + geom_line(aes(x=V4/1000, y=V1/60/60, colour=V2, group=V2), lwd=1.5) +
    geom_point(aes(x=V4/1000, y=V1/60/60, colour=V2), size=2.5) + 
    scale_colour_manual("Program", values=cols) +
    theme_bw() +
    scale_x_continuous("Number of cells (thousands)") +
    scale_y_continuous("Execution time (hours)") +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14))

ggsave(tplt, file="time_ncell.pdf", width=6, height=4)

dat <- dat[,c(8,10,11)]
dat2 <- dat2[,c(8,10,11)]
dat3$V10 <- dat3$V10 / 1000
dat3 <- dat3[,c(8,10,11)]
colnames(dat) <- c("prog", "x", "y")
colnames(dat2) <- c("prog", "x", "y")
colnames(dat3) <- c("prog", "x", "y")
dat$type <- "SNPs in panel"
dat2$type <- "Reads per cell"
dat3$type <- "Number of cells"
dat_all <- rbind(dat, dat2, dat3)

tdat <- tdat[,c(2,4,1)]
tdat2 <- tdat2[,c(2,4,1)]
tdat3 <- tdat3[,c(2,4,1)]
colnames(tdat) <- c("prog", "x", "y")
colnames(tdat2) <- c("prog", "x", "y")
colnames(tdat3) <- c("prog", "x", "y")
tdat2$x <- tdat2$x * 50
tdat$type <- "SNPs in panel"
tdat2$type <- "Reads per cell"
tdat3$type <- "Number of cells"
tdat3$x <- tdat3$x / 1000
tdat_all <- rbind(tdat, tdat2, tdat3)
tdat_all$y <- tdat_all$y/60/60

dat_all$cat <- "F score"
tdat_all$cat <- "Execution time (hours)"

all_all <- rbind(dat_all, tdat_all)

all_all <- all_all[order(all_all$prog),]

all_all$cat <- factor(all_all$cat, labels=c("F score", "Execution time (hours)"), levels=c("F score", "Execution time (hours)"))
all_all$type <- factor(all_all$type, labels=c("SNPs in panel (x1000)", "Reads per cell (x1000)", "Number of cells (x1000)"),
                       levels=c("SNPs in panel", "Reads per cell", "Number of cells"))

plt <- ggplot(all_all) + geom_line(aes(x=x, y=y, colour=prog, group=prog), lwd=1.5) + 
    geom_point(aes(x=x, y=y, colour=prog), size=3) + 
    theme_bw() + 
    facet_grid(cat~type, scales='free') + 
    scale_colour_manual("Program", values=cols) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          strip.background=element_blank(),
          strip.text.x=element_text(size=14),
          strip.text.y=element_text(size=14))

ggsave(plt, file="vcf_benchmark.pdf", width=9, height=5)


