#! /usr/bin/env Rscript

key1 <- read.table('combined_15_-1/key.txt', header=T)
key2a <- read.table('combined_20_-1/key.txt', header=T)
key2 <- read.table('combined_25_-1/key.txt', header=T)
key3 <- read.table('combined_35_-1/key.txt', header=T)
key4 <- read.table('combined_45_-1/key.txt', header=T)

key5 <- read.table('combined_20_1000000/key.txt', header=T)
key6 <- read.table('combined_20_5000000/key.txt', header=T)
key7 <- read.table('combined_20_10000000/key.txt', header=T)
key8 <- read.table('combined_20_20000000/key.txt', header=T)

map1 <- read.table('combined_15_-1/bcmap.txt')
map2a <- read.table('combined_20_-1/bcmap.txt')
map2 <- read.table('combined_25_-1/bcmap.txt')
map3 <- read.table('combined_35_-1/bcmap.txt')
map4 <- read.table('combined_45_-1/bcmap.txt')

map5 <- read.table('combined_20_1000000/bcmap.txt')
map6 <- read.table('combined_20_5000000/bcmap.txt')
map7 <- read.table('combined_20_10000000/bcmap.txt')
map8 <- read.table('combined_20_20000000/bcmap.txt')

key1$bc.orig <- map1$V1
key2a$bc.orig <- map2a$V1
key2$bc.orig <- map2$V1
key3$bc.orig <- map3$V1
key4$bc.orig <- map4$V1

key5$bc.orig <- map5$V1
key6$bc.orig <- map6$V1
key7$bc.orig <- map7$V1
key8$bc.orig <- map8$V1

ass15 <- read.table('combined_15_-1/species.assignments')
ass15 <- ass15[,c(1,2)]
colnames(ass15) <- c("bc", "k15")

f15 <- read.table('combined_15_-1/species.filt.assignments')
ass15$filt15 <- 0
ass15[which(ass15$bc %in% f15$V1),]$filt15 <- 1

ass15 <- merge(ass15, key1, by.x="bc", by.y="bc_uniq")

ass20 <- read.table('combined_20_-1/species.assignments')
ass20 <- ass20[,c(1,2)]
colnames(ass20) <- c("bc", "k20")

f20 <- read.table('combined_20_-1/species.filt.assignments')
ass20$filt20 <- 0
ass20[which(ass20$bc %in% f20$V1),]$filt20 <- 1

ass20 <- merge(ass20, key2a, by.x="bc", by.y="bc_uniq")

ass25 <- read.table('combined_25_-1/species.assignments')
ass25 <- ass25[,c(1,2)]
colnames(ass25) <- c("bc", "k25")

f25 <- read.table('combined_25_-1/species.filt.assignments')
ass25$filt25 <- 0
ass25[which(ass25$bc %in% f25$V1),]$filt25 <- 1

ass25 <- merge(ass25, key2, by.x="bc", by.y="bc_uniq")

ass35 <- read.table('combined_35_-1/species.assignments')
ass35 <- ass35[,c(1,2)]
colnames(ass35) <- c("bc", "k35")

f35 <- read.table('combined_35_-1/species.filt.assignments')
ass35$filt35 <- 0
ass35[which(ass35$bc %in% f35$V1),]$filt35 <- 1

ass35 <- merge(ass35, key3, by.x="bc", by.y="bc_uniq")

ass45 <- read.table('combined_45_-1/species.assignments')
ass45 <- ass45[,c(1,2)]
colnames(ass45) <- c("bc", "k45")

f45 <- read.table('combined_45_-1/species.filt.assignments')
ass45$filt45 <- 0
ass45[which(ass45$bc %in% f45$V1),]$filt45 <- 1

ass45 <- merge(ass45, key4, by.x="bc", by.y="bc_uniq")

ass20_1m <- read.table('combined_20_1000000/species.assignments')
ass20_1m <- ass20_1m[,c(1,2)]
colnames(ass20_1m) <- c("bc", "k20_1m")

f20_1m <- read.table('combined_20_1000000/species.filt.assignments')
ass20_1m$filt20_1m <- 0
ass20_1m[which(ass20_1m$bc %in% f20_1m$V1),]$filt20_1m <- 1

ass20_1m <- merge(ass20_1m, key5, by.x="bc", by.y="bc_uniq")

ass20_5m <- read.table('combined_20_5000000/species.assignments')
ass20_5m <- ass20_5m[,c(1,2)]
colnames(ass20_5m) <- c("bc", "k20_5m")

f20_5m <- read.table('combined_20_5000000/species.filt.assignments')
ass20_5m$filt20_5m <- 0
ass20_5m[which(ass20_5m$bc %in% f20_5m$V1),]$filt20_5m <- 1

ass20_5m <- merge(ass20_5m, key6, by.x="bc", by.y="bc_uniq")

ass20_10m <- read.table('combined_20_10000000/species.assignments')
ass20_10m <- ass20_10m[,c(1,2)]
colnames(ass20_10m) <- c("bc", "k20_10m")

f20_10m <- read.table('combined_20_10000000/species.filt.assignments')
ass20_10m$filt20_10m <- 0
ass20_10m[which(ass20_10m$bc %in% f20_10m$V1),]$filt20_10m <- 1

ass20_10m <- merge(ass20_10m, key7, by.x="bc", by.y="bc_uniq")

ass20_20m <- read.table('combined_20_20000000/species.assignments')
ass20_20m <- ass20_20m[,c(1,2)]
colnames(ass20_20m) <- c("bc", "k20_20m")

f20_20m <- read.table('combined_20_20000000/species.filt.assignments')
ass20_20m$filt20_20m <- 0
ass20_20m[which(ass20_20m$bc %in% f20_20m$V1),]$filt20_20m <- 1

ass20_20m <- merge(ass20_20m, key8, by.x="bc", by.y="bc_uniq")

merged <- merge(ass15[,-c(1)], ass20[,-c(1)], by=c("true", "bc.orig", "filtered_set"), all.x=TRUE, all.y=TRUE)
merged <- merge(merged, ass25[,-c(1)], by=c("true", "bc.orig", "filtered_set"), all.x=TRUE, all.y=TRUE)
merged <- merge(merged, ass35[,-c(1)], by=c("true", "bc.orig", "filtered_set"), all.x=TRUE, all.y=TRUE)
merged <- merge(merged, ass45[,-c(1)], by=c("true", "bc.orig", "filtered_set"), all.x=TRUE, all.y=TRUE)
merged <- merge(merged, ass20_1m[,-c(1)], by=c("true", "bc.orig", "filtered_set"), all.x=TRUE, all.y=TRUE)
merged <- merge(merged, ass20_5m[,-c(1)], by=c("true", "bc.orig", "filtered_set"), all.x=TRUE, all.y=TRUE)
merged <- merge(merged, ass20_10m[,-c(1)], by=c("true", "bc.orig", "filtered_set"), all.x=TRUE, all.y=TRUE)
merged <- merge(merged, ass20_20m[,-c(1)], by=c("true", "bc.orig", "filtered_set"), all.x=TRUE, all.y=TRUE)

##### Table heading explanation #####
# sens: sensitivity; the fraction of "true" calls correctly assigned by CellBouncer. This includes
#    empty droplets, as these barcodes were neither filtered by CellBouncer nor the publication.
#    normalized by the number of cells with a true identity (all cells)
# sens_filt: fraction of publication-filtered assignments correctly assigned by CellBouncer
# spec: specificity; the fraction of all CellBouncer assignments that match the true (publication)
#    assignment (NOTE: this will be low because it includes data from empty/low quality droplets!)
#    normalized by the number of cells with a CellBouncer-assigned identity
# spec_filt: fraction of CellBouncer-filtered assignments that match truth from publication
#    (this measures how good CellBouncer was at separating true cells from empty droplets)
# npresent: Total number of barcodes assigned an identity by CellBouncer (includes empty droplets)
# nfilt: Total number of CellBouncer-filtered barcodes (excludes CellBouncer's guess of empty droplets)
# npub: Total number of publication-filtered barcodes assigned an identity
# f: (2*sens*spec)/(sens + spec)
# f_filt: (2*sens_filt*spec_filt)/(sens_filt + spec_filt)

k15_sens <- length(rownames(merged[which(merged$k15==merged$true),]))/length(rownames(merged[which(!is.na(merged$true)),]))
k15_sens_filt <- length(rownames(merged[which(merged$k15==merged$true & merged$filtered_set==1),])) / 
    length(rownames(merged[which(merged$filtered_set==1),]))
k15_spec <- length(rownames(merged[which(merged$k15==merged$true),]))/length(rownames(merged[which(!is.na(merged$k15)),]))
k15_spec_filt <- length(rownames(merged[which(merged$k15==merged$true & merged$filt15==1),])) /
    length(rownames(merged[which(merged$filt15==1),]))
k15_npresent <- length(rownames(merged[which(!is.na(merged$k15)),]))
k15_nfilt <- length(rownames(merged[which(merged$filt15==1),]))
k15_npub <- length(rownames(merged[which(merged$filtered_set==1 & !is.na(merged$k15)),]))
k15_f <- 2*(k15_sens*k15_spec)/(k15_sens+k15_spec)
k15_f_filt <- 2*(k15_sens_filt*k15_spec_filt)/(k15_sens_filt+k15_spec_filt)

k20_sens <- length(rownames(merged[which(merged$k20==merged$true),]))/length(rownames(merged[which(!is.na(merged$true)),]))
k20_sens_filt <- length(rownames(merged[which(merged$k20==merged$true & merged$filtered_set==1),])) / 
    length(rownames(merged[which(merged$filtered_set==1),]))
k20_spec <- length(rownames(merged[which(merged$k20==merged$true),]))/length(rownames(merged[which(!is.na(merged$k20)),]))
k20_spec_filt <- length(rownames(merged[which(merged$k20==merged$true & merged$filt20==1),])) /
    length(rownames(merged[which(merged$filt20==1),]))
k20_npresent <- length(rownames(merged[which(!is.na(merged$k20)),]))
k20_nfilt <- length(rownames(merged[which(merged$filt20==1),]))
k20_npub <- length(rownames(merged[which(merged$filtered_set==1 & !is.na(merged$k20)),]))
k20_f <- 2*(k20_sens*k20_spec)/(k20_sens+k20_spec)
k20_f_filt <- 2*(k20_sens_filt*k20_spec_filt)/(k20_sens_filt+k20_spec_filt)

k25_sens <- length(rownames(merged[which(merged$k25==merged$true),]))/length(rownames(merged[which(!is.na(merged$true)),]))
k25_sens_filt <- length(rownames(merged[which(merged$k25==merged$true & merged$filtered_set==1),])) / 
    length(rownames(merged[which(merged$filtered_set==1),]))
k25_spec <- length(rownames(merged[which(merged$k25==merged$true),]))/length(rownames(merged[which(!is.na(merged$k25)),]))
k25_spec_filt <- length(rownames(merged[which(merged$k25==merged$true & merged$filt25==1),])) /
    length(rownames(merged[which(merged$filt25==1),]))
k25_npresent <- length(rownames(merged[which(!is.na(merged$k25)),]))
k25_nfilt <- length(rownames(merged[which(merged$filt25==1),]))
k25_npub <- length(rownames(merged[which(merged$filtered_set==1 & !is.na(merged$k25)),]))
k25_f <- 2*(k25_sens*k25_spec)/(k25_sens+k25_spec)
k25_f_filt <- 2*(k25_sens_filt*k25_spec_filt)/(k25_sens_filt+k25_spec_filt)

k35_sens <- length(rownames(merged[which(merged$k35==merged$true),]))/length(rownames(merged[which(!is.na(merged$true)),]))
k35_sens_filt <- length(rownames(merged[which(merged$k35==merged$true & merged$filtered_set==1),])) / 
    length(rownames(merged[which(merged$filtered_set==1),]))
k35_spec <- length(rownames(merged[which(merged$k35==merged$true),]))/length(rownames(merged[which(!is.na(merged$k35)),]))
k35_spec_filt <- length(rownames(merged[which(merged$k35==merged$true & merged$filt35==1),])) /
    length(rownames(merged[which(merged$filt35==1),]))
k35_npresent <- length(rownames(merged[which(!is.na(merged$k35)),]))
k35_nfilt <- length(rownames(merged[which(merged$filt35==1),]))
k35_npub <- length(rownames(merged[which(merged$filtered_set==1 & !is.na(merged$k35)),]))
k35_f <- 2*(k35_sens*k35_spec)/(k35_sens+k35_spec)
k35_f_filt <- 2*(k35_sens_filt*k35_spec_filt)/(k35_sens_filt+k35_spec_filt)

k45_sens <- length(rownames(merged[which(merged$k45==merged$true),]))/length(rownames(merged[which(!is.na(merged$true)),]))
k45_sens_filt <- length(rownames(merged[which(merged$k45==merged$true & merged$filtered_set==1),])) / 
    length(rownames(merged[which(merged$filtered_set==1),]))
k45_spec <- length(rownames(merged[which(merged$k45==merged$true),]))/length(rownames(merged[which(!is.na(merged$k45)),]))
k45_spec_filt <- length(rownames(merged[which(merged$k45==merged$true & merged$filt45==1),])) /
    length(rownames(merged[which(merged$filt45==1),]))
k45_npresent <- length(rownames(merged[which(!is.na(merged$k45)),]))
k45_nfilt <- length(rownames(merged[which(merged$filt45==1),]))
k45_npub <- length(rownames(merged[which(merged$filtered_set==1 & !is.na(merged$k45)),]))
k45_f <- 2*(k45_sens*k45_spec)/(k45_sens+k45_spec)
k45_f_filt <- 2*(k45_sens_filt*k45_spec_filt)/(k45_sens_filt+k45_spec_filt)

k20_1m_sens <- length(rownames(merged[which(merged$k20_1m==merged$true),]))/length(rownames(merged[which(!is.na(merged$true)),]))
k20_1m_sens_filt <- length(rownames(merged[which(merged$k20_1m==merged$true & merged$filtered_set==1),])) / 
    length(rownames(merged[which(merged$filtered_set==1),]))
k20_1m_spec <- length(rownames(merged[which(merged$k20_1m==merged$true),]))/length(rownames(merged[which(!is.na(merged$k20_1m)),]))
k20_1m_spec_filt <- length(rownames(merged[which(merged$k20_1m==merged$true & merged$filt20_1m==1),])) /
    length(rownames(merged[which(merged$filt20_1m==1),]))
k20_1m_npresent <- length(rownames(merged[which(!is.na(merged$k20_1m)),]))
k20_1m_nfilt <- length(rownames(merged[which(merged$filt20_1m==1),]))
k20_1m_npub <- length(rownames(merged[which(merged$filtered_set==1 & !is.na(merged$k20_1m)),]))
k20_1m_f <- 2*(k20_1m_sens*k20_1m_spec)/(k20_1m_sens+k20_1m_spec)
k20_1m_f_filt <- 2*(k20_1m_sens_filt*k20_1m_spec_filt)/(k20_1m_sens_filt+k20_1m_spec_filt)

k20_5m_sens <- length(rownames(merged[which(merged$k20_5m==merged$true),]))/length(rownames(merged[which(!is.na(merged$true)),]))
k20_5m_sens_filt <- length(rownames(merged[which(merged$k20_5m==merged$true & merged$filtered_set==1),])) / 
    length(rownames(merged[which(merged$filtered_set==1),]))
k20_5m_spec <- length(rownames(merged[which(merged$k20_5m==merged$true),]))/length(rownames(merged[which(!is.na(merged$k20_5m)),]))
k20_5m_spec_filt <- length(rownames(merged[which(merged$k20_5m==merged$true & merged$filt20_5m==1),])) /
    length(rownames(merged[which(merged$filt20_5m==1),]))
k20_5m_npresent <- length(rownames(merged[which(!is.na(merged$k20_5m)),]))
k20_5m_nfilt <- length(rownames(merged[which(merged$filt20_5m==1),]))
k20_5m_npub <- length(rownames(merged[which(merged$filtered_set==1 & !is.na(merged$k20_5m)),]))
k20_5m_f <- 2*(k20_5m_sens*k20_5m_spec)/(k20_5m_sens+k20_5m_spec)
k20_5m_f_filt <- 2*(k20_5m_sens_filt*k20_5m_spec_filt)/(k20_5m_sens_filt+k20_5m_spec_filt)

k20_10m_sens <- length(rownames(merged[which(merged$k20_10m==merged$true),]))/length(rownames(merged[which(!is.na(merged$true)),]))
k20_10m_sens_filt <- length(rownames(merged[which(merged$k20_10m==merged$true & merged$filtered_set==1),])) / 
    length(rownames(merged[which(merged$filtered_set==1),]))
k20_10m_spec <- length(rownames(merged[which(merged$k20_10m==merged$true),]))/length(rownames(merged[which(!is.na(merged$k20_10m)),]))
k20_10m_spec_filt <- length(rownames(merged[which(merged$k20_10m==merged$true & merged$filt20_10m==1),])) /
    length(rownames(merged[which(merged$filt20_10m==1),]))
k20_10m_npresent <- length(rownames(merged[which(!is.na(merged$k20_10m)),]))
k20_10m_nfilt <- length(rownames(merged[which(merged$filt20_10m==1),]))
k20_10m_npub <- length(rownames(merged[which(merged$filtered_set==1 & !is.na(merged$k20_10m)),]))
k20_10m_f <- 2*(k20_10m_sens*k20_10m_spec)/(k20_10m_sens+k20_10m_spec)
k20_10m_f_filt <- 2*(k20_10m_sens_filt*k20_10m_spec_filt)/(k20_10m_sens_filt+k20_10m_spec_filt)

k20_20m_sens <- length(rownames(merged[which(merged$k20_20m==merged$true),]))/length(rownames(merged[which(!is.na(merged$true)),]))
k20_20m_sens_filt <- length(rownames(merged[which(merged$k20_20m==merged$true & merged$filtered_set==1),])) / 
    length(rownames(merged[which(merged$filtered_set==1),]))
k20_20m_spec <- length(rownames(merged[which(merged$k20_20m==merged$true),]))/length(rownames(merged[which(!is.na(merged$k20_20m)),]))
k20_20m_spec_filt <- length(rownames(merged[which(merged$k20_20m==merged$true & merged$filt20_20m==1),])) /
    length(rownames(merged[which(merged$filt20_20m==1),]))
k20_20m_npresent <- length(rownames(merged[which(!is.na(merged$k20_20m)),]))
k20_20m_nfilt <- length(rownames(merged[which(merged$filt20_20m==1),]))
k20_20m_npub <- length(rownames(merged[which(merged$filtered_set==1 & !is.na(merged$k20_20m)),]))
k20_20m_f <- 2*(k20_20m_sens*k20_20m_spec)/(k20_20m_sens+k20_20m_spec)
k20_20m_f_filt <- 2*(k20_20m_sens_filt*k20_20m_spec_filt)/(k20_20m_sens_filt+k20_20m_spec_filt)

df <- data.frame(k=c(15,20,25,35,45,20,20,20,20),
    samp=c("all", "all", "all", "all", "all", "1m", "5m", "10m", "20m"),
    f=c(k15_f, k20_f, k25_f, k35_f, k45_f, k20_1m_f, k20_5m_f, k20_10m_f, k20_20m_f),
    f_filt=c(k15_f_filt, k20_f_filt, k25_f_filt, k35_f_filt, k45_f_filt, k20_1m_f_filt, k20_5m_f_filt, k20_10m_f_filt, k20_20m_f_filt),
    sens=c(k15_sens, k20_sens, k25_sens, k35_sens, k45_sens, k20_1m_sens, k20_5m_sens, k20_10m_sens, k20_20m_sens),
    sens_filt=c(k15_sens_filt, k20_sens_filt, k25_sens_filt, k35_sens_filt, k45_sens_filt, k20_1m_sens_filt, k20_5m_sens_filt, k20_10m_sens_filt, k20_20m_sens_filt),
    spec=c(k15_spec, k20_spec, k25_spec, k35_spec, k45_spec, k20_1m_spec, k20_5m_spec, k20_10m_spec, k20_20m_spec),
    spec_filt=c(k15_spec_filt, k20_spec_filt, k25_spec_filt, k35_spec_filt, k45_spec_filt, k20_1m_spec_filt, k20_5m_spec_filt, k20_10m_spec_filt, k20_20m_spec_filt),
    npresent=c(k15_npresent, k20_npresent, k25_npresent, k35_npresent, k45_npresent, k20_1m_npresent, k20_5m_npresent, k20_10m_npresent, k20_20m_npresent),
    nfilt=c(k15_nfilt, k20_nfilt, k25_nfilt, k35_nfilt, k45_nfilt, k20_1m_nfilt, k20_5m_nfilt, k20_10m_nfilt, k20_20m_nfilt),
    npub=c(k15_npub, k20_npub, k25_npub, k35_npub, k45_npub, k20_1m_npub, k20_5m_npub, k20_10m_npub, k20_20m_npub))

write.table(df, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

write.table(df, file="kmercomp_summary.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(merged, file="kmercomp_data.txt", sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

