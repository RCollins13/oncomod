#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Compute median value per sample from a two-column tsv of sample, value


# Single command line argument: input.tsv
tsv.in <- as.character(commandArgs(trailingOnly=T)[1])

# Load input data and coerce second column to numeric
data <- read.table(tsv.in, header=F, sep="\t")
colnames(data) <- c("sample", "value")
data$value <- as.numeric(data$value)

# Compute median value per sample
samples <- sort(unique(data$sample))
medians <- sapply(samples, function(sample){
  median(data[which(data$sample == sample), "value"], na.rm=T)
})
res <- data.frame("sample" = samples, "median" = medians)

# Write results to stdout
write.table(res, col.names=F, row.names=F, sep="\t", quote=F)
