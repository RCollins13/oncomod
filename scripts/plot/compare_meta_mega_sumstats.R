#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Compare P-values and effect sizes between meta-analysis and pooled mega-analysis


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(rCNV2, quietly=TRUE)
require(OncoModR, quietly=TRUE)
OncoModR::load.constants("all")


##################
# Data Functions #
##################
load.sumstats <- function(tsv.in, pval.column="p", beta.column="beta"){
  ss <- read.table(tsv.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(ss)[1] <- gsub("#", "", colnames(ss)[1])
  ss[, pval.column] <- -log10(ss[, pval.column])
  ss <- ss[, c("somatic", "germline", beta.column, pval.column)]
  ss[complete.cases(ss), ]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Compare summary statistics between ",
                                           "meta-analysis and pooled mega-",
                                           "analysis models", sep=""))
parser$add_argument("--meta-sumstats", metavar=".tsv", type="character",
                    help="meta-analysis summary statistics .tsv", required=TRUE)
parser$add_argument("--pooled-sumstats", metavar=".tsv", type="character",
                    help="pooled mega-analysis summary statistics .tsv", required=TRUE)
parser$add_argument("--cancer", metafar="string", type="character",
                    help="(Optional) Cancer type abbreviation to color plots by cancer")
parser$add_argument("--out-prefix", metavar="path|string", type="character",
                    default="sumstat_comparison", help="file prefix for output plots")
args <- parser$parse_args()

# # DEV:
# args <- list("meta_sumstats" = "~/scratch/CRAD.meta.sumstats.tsv.gz",
#              "pooled_sumstats" = "~/scratch/pooled.CRAD.sumstats.tsv.gz",
#              "cancer" = "CRAD",
#              "out_prefix" = "~/scratch/sumstat_comparison_test")

# Load P-values and odds ratios from meta & pooled summary statistics
meta <- load.sumstats(args$meta_sumstats)
pooled <- load.sumstats(args$pooled_sumstats)

# Merge meta & pooled tests, retaining only tests found in both
data <- merge(meta, pooled, by=c("somatic", "germline"), all=F,
              suffixes=c(".meta", ".pooled"), sort=F)
data <- data[complete.cases(data), ]

# Set cancer-type specific palette, if optioned
if(!is.null(args$cancer)){
  pal <- colorRampPalette(cancer.palettes[[toupper(args$cancer)]][3:7])
}else{
  pal <- NULL
}

# Plot P-values
max.p <- ceiling(max(data[, paste("p", c("meta", "pooled"), sep=".")], na.rm=T))
png(paste(args$out_prefix, "meta_vs_pooled.p_values.png", sep="."),
    height=3.25*300, width=3.25*300, res=300)
par(mar=c(2.75, 2.75, 0.5, 0.5), bty="n")
plot(NA, xlim=c(0, max.p), ylim=c(0, max.p), xaxt="n", yaxt="n",
     xlab="", ylab="", xaxs="i", yaxs="i")
abline(0, 1, lty=5)
clean.axis(1, infinite=TRUE, title=bquote("Pooled" ~ -log[10](italic(P))))
clean.axis(2, infinite=TRUE, title=bquote("Meta-Analysis" ~ -log[10](italic(P))))
rCNV2::dens.scatter(x=data$p.pooled, y=data$p.meta, palette=pal,
                    pt.cex=0.2, plot.cor=F, add=T)
dev.off()

# Plot effect sizes
max.beta <- min(c(5, ceiling(max(abs(data[, paste("beta", c("meta", "pooled"), sep=".")]), na.rm=T))))
png(paste(args$out_prefix, "meta_vs_pooled.betas.png", sep="."),
    height=3.25*300, width=3.25*300, res=300)
par(mar=c(2.75, 2.75, 0.5, 0.5), bty="n")
plot(NA, xlim=c(-max.beta, max.beta), ylim=c(-max.beta, max.beta), xaxt="n", yaxt="n",
     xlab="", ylab="", xaxs="i", yaxs="i")
abline(0, 1, lty=5)
clean.axis(1, infinite=TRUE, title=bquote("Pooled" ~ ln("OR")))
clean.axis(2, infinite=TRUE, title=bquote("Meta-Analysis" ~ ln("OR")))
rCNV2::dens.scatter(x=data$beta.pooled, y=data$beta.meta, palette=pal,
                    pt.cex=0.2, plot.cor=F, add=T)
dev.off()
