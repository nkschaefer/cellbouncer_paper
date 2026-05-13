#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1){
    write("ERROR: please provide output.ps file, index", stderr())
    q()    
}

# Read lines from dot.ps
lines <- readLines(args[1])

# Extract lines with 'ubox'
ubox_lines <- grep("ubox", lines, value = TRUE)
ubox_lines <- ubox_lines[grep("%", ubox_lines, fixed=T, invert=T)] 
ubox_lines <- ubox_lines[grep("{", ubox_lines, fixed=T, invert=T)]

# Parse the data
data <- do.call(rbind, strsplit(ubox_lines, " "))
data <- data.frame(
  i = as.integer(data[, 1]),
  j = as.integer(data[, 2]),
  sqrt_p = as.numeric(data[, 3])
)

# Add real probability
data$p <- data$sqrt_p^2

allp <- rbind(data.frame(x=data$i, p=data$p), data.frame(x=data$j, p=data$p))

idx <- as.numeric(args[2])
write.table(max(allp[which(allp$x==idx),]$p), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

