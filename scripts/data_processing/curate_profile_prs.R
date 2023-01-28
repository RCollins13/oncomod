#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Curate selected PRS for PROFILE samples

# Set options
options(stringsAsFactors=F, scipen=1000)

# Read command-line arguments
args <- commandArgs(trailingOnly=TRUE)

# Load PRS matrix
prs <- read.table(args[1], sep=" ", check.names=F, header=T, comment.char="")

# Subset to list of samples included in analysis
samples <- read.table(args[2], header=F, sep="")[, 1]
prs <- prs[which(prs[, 1] %in% samples), ]

# Subset to PRS of interest
keepers <- read.table(args[3], header=F, sep="\t")[, 2]
prs <- prs[, c(1, which(colnames(prs) %in% keepers))]

# Transpose and write to outfile
rownames(prs) <- prs[, 1]
prs[, 1] <- NULL
VID <- colnames(prs)
prs <- t(prs)
prs <- as.data.frame(cbind(VID, prs))
pres <- prs[, c("VID", setdiff(colnames(prs), "VID"))]
colnames(prs)[1] <- paste("#", colnames(prs)[1], sep="")
write.table(prs, args[4], col.names=T, row.names=F, sep="\t", quote=F)
