#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Visualize FGFR4 haplotypes


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(OncoModR, quietly=TRUE)
require(viridis, quietly=TRUE)
load.constants("all")


##################
# Data Functions #
##################
# Load an LD matrix
load.ld <- function(ld.in){
  # Read LD data
  ld.long <- read.table(ld.in, sep="\t", header=T)
  vids <- sort(unique(unlist(ld.long[, 1:2])))

  # Transform to square all-by-all matrix
  ld.matrix <- as.data.frame(do.call("rbind", lapply(vids, function(vid.r){
    unlist(sapply(vids, function(vid.c){
      if(vid.c == vid.r){return(1)}
      hit.counts <- apply(ld.long[, 1:2], 1, function(rv){length(intersect(rv, c(vid.r, vid.c)))})
      hit <- which(hit.counts == 2)
      if(length(hit) == 0){0}else{ld.long[hit, 3]}
    }))
  })))
  colnames(ld.matrix) <- rownames(ld.matrix) <- vids

  # Cluster rows and columns
  order <- hclust(dist(ld.matrix))$order
  ld.matrix[order, order]
}


######################
# Plotting Functions #
######################
# Plot LD heatmap
plot.ld.matrix <- function(ld.matrix, pal.fx=viridis, label.vars=NULL,
                           label.colors=NULL){
  # Get plot dimensions
  n.vars <- nrow(ld.matrix)
  plot.dat <- floor(ld.matrix * 100) + 1
  pal <- pal.fx(101)

  # Prep plot area
  par(mar=c(0.2, 2, 2, 0.2), bty="n")
  plot(NA, xlim=c(0, n.vars), ylim=c(n.vars, 0),
       xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="")

  # Color cells
  sapply(1:n.vars, function(r){
    rect(xleft=(1:n.vars)-1, xright=1:n.vars, ybottom=r-1, ytop=r,
         col=pal[as.numeric(plot.dat[r, ])], lwd=0.2, xpd=T)
  })

  # Label selected variants, if optioned
  if(!is.null(label.vars)){
    if(is.null(label.colors)){
      label.colors <- rep("black", length(label.vars))
      names(label.colors) <- names(label.vars)
    }
    label.hits <- intersect(names(label.vars), rownames(ld.matrix))
    label.idxs <- sapply(label.hits, function(vid){})
    sapply(label.hits, function(vid){
      idx <- which(rownames(ld.matrix) == vid)
      axis(2, idx-0.5, las=2, cex.axis=3/6, col.axis=label.colors[vid],
           labels=label.vars[vid], line=-0.9, tick=F)
      axis(3, idx-0.5, las=2, cex.axis=3/6, col.axis=label.colors[vid],
           labels=label.vars[vid], line=-0.9, tick=F)
    })
  }
}


###########
# RScript #
###########
# Four positional arguments
args <- commandArgs(trailingOnly=TRUE)
coding.ld.in <- as.character(args[1])
common.ld.in <- as.character(args[2])
variant.sets.in <- as.character(args[3])
dir.out <- as.character(args[4])

# # DEV
# coding.ld.in <- "~/scratch/all_cohorts.FGFR4.CRAD.coding.ld.tsv.gz"
# common.ld.in <- "~/scratch/all_cohorts.FGFR4.CRAD.common_biallelic.ld.tsv.gz"
# variant.sets.in <- "~/scratch/FGFR4.coding_variants.sets.tsv"
# dir.out <- "~/scratch/"

# Load variant sets
vinfo <- read.table(variant.sets.in, header=F, sep="\t")
vinfo$key <- sapply(vinfo[, 2], function(ss){
  unique(sapply(unlist(strsplit(ss, split=",")), function(s){
    sparts <- unlist(strsplit(s, split="_"))
    paste(sparts[(length(sparts)-3):length(sparts)], collapse="_")
  }))
})
vinfo$csq <- "missense"
vinfo$csq[grep("=$", vinfo[, 1])] <- "synonymous"
vinfo$csq[grep("Ter$|Ter[0-9]+$", vinfo[, 1])] <- "lof"

# Load LD
coding.ld <- load.ld(coding.ld.in)
keep.coding.idxs <- which(rownames(coding.ld) %in% unlist(vinfo$key))
coding.ld <- coding.ld[keep.coding.idxs, keep.coding.idxs]
common.ld <- load.ld(common.ld.in)

# Custom labels
label.sets <- c("ENST00000292408_p.Gly388Arg" = "G388R",
                "ENST00000292408_p.Pro136Leu" = "P136L",
                "ENST00000292408_p.Arg54=" = "R54R",
                "ENST00000292408_p.Ser137=" = "S137S",
                "ENST00000292408_p.Asp575Asn" = "D575N",
                "ENST00000292408_p.Val10Ile" = "V10I",
                "ENST00000292408_p.Arg234=" = "R234R")
label.set.idxs <- sapply(names(label.sets), function(sid){which(vinfo[, 1] == sid)})
label.vars <- label.sets
label.colors <- csq.colors[vinfo$csq[label.set.idxs]]
names(label.colors) <- names(label.vars) <- unlist(vinfo[label.set.idxs, "key"])

# Plot coding LD
pdf(paste(dir.out, "FGFR4.LD.coding.pdf", sep="/"),
    height=4, width=4)
plot.ld.matrix(coding.ld, label.vars=label.vars, label.colors=label.colors)
dev.off()

# Plot common LD
pdf(paste(dir.out, "FGFR4.LD.common.pdf", sep="/"),
    height=4, width=4)
plot.ld.matrix(common.ld, label.vars=label.vars, label.colors=label.colors)
dev.off()

