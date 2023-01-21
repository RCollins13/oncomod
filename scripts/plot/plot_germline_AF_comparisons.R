#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Plot comparisons of germline allele frequencies between all cohorts


#########
# Setup #
#########
# Load necessary libraries and constants
require(RASMod, quietly=TRUE)
require(rCNV2, quietly=TRUE)
require(argparse, quietly=TRUE)
RASMod::load.constants("all")


##################
# Data functions #
##################
# Load vcf2bed output for a single cohort
load.vcf2bed <- function(file, suffix=NULL){
  df <- read.table(file, check.names=F, comment.char="", sep="\t", header=T)
  colnames(df)[1] <- gsub("^#", "", colnames(df)[1])
  suffix.idxs <- grep("_AF$", colnames(df))
  keep.col.idxs <- sort(unique(c(1:3, 5,
                                 grep("^gnomAD_", colnames(df)),
                                 suffix.idxs)))
  if(!is.null(suffix)){
    colnames(df)[suffix.idxs] <- paste(colnames(df)[suffix.idxs], suffix, sep=".")
  }
  df[, keep.col.idxs]
}


######################
# Plotting functions #
######################
# Plot a density-colored scatterplot of AF correlations between two cohorts
plot.AF.scatter <- function(x, y, pair.names, pop){
  prep.plot.area(xlims=c(0, 1), ylims=c(0, 1), parmar=c(2.5, 2.5, 1.25, 1.25))
  dens.scatter(x=x, y=y, add=T, pt.cex=0.3,
               palette=colorRampPalette(get(paste(pop, "colors", sep="."))))
  sapply(1:2, function(k){
    clean.axis(k, title=paste(pair.names[k], "AF"), infinite=T)
  })
  mtext(paste(pair.names[1], " vs. ", pair.names[2], " (",
              pop.names.short[pop], ")", sep=""), 3, line=0.25)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Summarize cohort and patient phenotype metadata")
parser$add_argument("--bed", metavar="BED", type="character", action="append",
                    help="single-cohort vcf2bed output", required=TRUE)
parser$add_argument("--name", metavar="string", type="character", action="append",
                    help=paste("names for each --bed input, provided",
                               "in the same order as --bed"))
parser$add_argument('--out-prefix', metavar='path', type="character", nargs=1,
                    help='file prefix for all plots')
args <- parser$parse_args()

# # DEV
# args <- list("bed"=c("~/scratch/PROFILE.germline_variants.bed.gz",
#                        "~/scratch/TCGA.germline_variants.bed.gz"),
#              "name"=c("PROFILE", "TCGA"),
#              "out_prefix"="~/scratch/AF_cor_test")

# Sanity-check lengths of --bed and --name
if(is.null(args$name)){
  args$name <- paste("cohort", 1:length(args$bed), sep="")
}else{
  if(length(args$bed) != length(args$name)){
    stop("Different number of --bed and --name arguments provided. Exiting.")
  }
}
n.cohorts <- length(args$bed)

# Load data for each cohort
dat.list <- apply(data.frame(args$bed, args$name), 1,
                    function(rvals){
                      load.vcf2bed(rvals[1], suffix=rvals[2])
                    })
names(dat.list) <- args$name

# Merge data by overlapping columns
dat <- Reduce(function(x, y){merge(x, y, all=T, sort=F)}, dat.list)


############################
# Inter-cohort comparisons #
############################
apply(combn(args$name, 2), 2, function(pair.names){
  # One plot per ancestry
  sapply(names(pop.names.short), function(pop){
    col.names <- paste(pop, "_AF.", pair.names, sep="")
    if(!all(col.names %in% colnames(dat))){return()}
    sub.dat <- dat[, col.names]
    sub.dat <- sub.dat[complete.cases(sub.dat), ]
    if(nrow(sub.dat) == 0){return()}
    png(paste(args$out_prefix, "intercohort_AF_comparison",
              pair.names[1], pair.names[2], pop, "png", sep="."),
        height=3*300, width=3*300, res=300)

    dev.off()
  })
})




