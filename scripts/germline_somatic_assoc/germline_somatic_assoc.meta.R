#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Conduct germline-somatic association meta-analyses of multiple cohorts


#########
# Setup #
#########
# Load necessary libraries and constants
require(RASMod, quietly=TRUE)
require(argparse, quietly=TRUE)
RASMod::load.constants("names")


##################
# Data functions #
##################
# Load association statistics from a single cohort
load.stats.single <- function(file, suffix=NULL){
  df <- read.table(file, check.names=F, comment.char="", sep="\t", header=T)
  colnames(df)[1] <- gsub("^#", "", colnames(df)[1])
  if(!is.null(suffix)){
    colnames(df)[-c(1:2)] <- paste(colnames(df)[-c(1:2)], suffix, sep=".")
  }
  df[, -c(1:2)] <- apply(df[, -c(1:2)], 2, as.numeric)
  return(df)
}

# Merge a list of single-cohort association statistics, compute weighted frequencies,
# and run inverse-variance weighted meta-analysis of effect sizes
merge.stats <- function(stats.list){
  # Iteratively outer join all stats
  merged <- Reduce(function(x, y){merge(x, y, all=T, by=c("somatic", "germline"), sort=F)},
                   stats.list)

  # Gather information on non-NA cohorts
  cohort.info <- as.data.frame(t(apply(merged[, grep("^z\\.", colnames(merged))], 1,
                                       function(zscores){
                                         idxs <- which(!is.na(zscores) & !is.infinite(zscores))
                                         cohorts <- names(stats.list)[idxs]
                                         c(length(cohorts), paste(sort(cohorts), collapse=","))
                                       })))
  colnames(cohort.info) <- c("N_cohorts", "cohorts")

  # Compute weighted frequencies
  wfreqs <- as.data.frame(t(apply(merged, 1, function(rvals){
    samps <- as.numeric(rvals[grep("^samples\\.", names(rvals))])
    weights <- sqrt(samps)
    som.ac <- as.numeric(rvals[grep("^somatic_AC\\.", names(rvals))])
    som.freqs <- som.ac / samps
    som.freq.w <- weighted.mean(som.freqs, weights, na.rm=T)
    yes_som.germ.ac <- as.numeric(rvals[grep("^yes_somatic.germline_AC\\.", names(rvals))])
    yes_som.germ.freqs <- yes_som.germ.ac / (2 * samps)
    yes_som.germ.freq.w <- weighted.mean(yes_som.germ.freqs, weights, na.rm=T)
    no_som.germ.ac <- as.numeric(rvals[grep("^no_somatic.germline_AC\\.", names(rvals))])
    no_som.germ.freqs <- no_som.germ.ac / (2 * samps)
    no_som.germ.freq.w <- weighted.mean(no_som.germ.freqs, weights, na.rm=T)
    return(c(som.freq.w, yes_som.germ.freq.w, no_som.germ.freq.w))
  })))
  colnames(wfreqs) <- c("somatic_freq", "yes_somatic.germline_AF", "no_somatic.germline_AF")

  # Add cohort info & weighted frequencies to merged dataframe and drop all single-cohort freq info
  cols.to.drop <- grep("^samples\\.|^somatic_AC\\.|\\.germline_AC\\.", colnames(merged))
  merged <- cbind(merged[, -cols.to.drop], cohort.info, wfreqs)

  # Meta-analyze effect sizes
  meta.stats <- as.data.frame(t(apply(merged, 1, function(rvals){
    betas <- as.numeric(rvals[grep("^beta\\.", names(rvals))])
    ses <- as.numeric(rvals[grep("^beta_SE\\.", names(rvals))])
    ivw.meta(betas, ses)
  })))
  colnames(meta.stats) <- c("beta", "beta_SE")
  meta.stats$z <- meta.stats$beta / meta.stats$beta_SE
  meta.stats$p <- 2*pnorm(-abs(meta.stats$z))

  # Add meta-analysis stats to merged results and drop all single-cohort stats
  cols.to.drop <- grep("^beta\\.|^beta_SE\\.|^z\\.|^p\\.", colnames(merged))
  merged <- cbind(merged[, -cols.to.drop], meta.stats)

  # Sort rows and reorder columns
  merged[with(merged, order(somatic, germline)),
         c("somatic", "germline", "somatic_freq", "yes_somatic.germline_AF",
           "no_somatic.germline_AF", "N_cohorts", "cohorts",
           "beta", "beta_SE", "z", "p")]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Conduct multi-cohort germline-somatic",
                                           "association meta-analysis"))
parser$add_argument("--stats", metavar=".tsv", type="character", action="append",
                    help="single-cohort association stats", required=TRUE)
parser$add_argument("--name", metavar="string", type="character", action="append",
                    help=paste("names for each --stats input, provided",
                               "in the same order as --stats"))
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="output .tsv file for meta-analysis statistics")
args <- parser$parse_args()

# # DEV
# args <- list("stats"=c("~/scratch/PROFILE.SKCM.sumstats.tsv.gz",
#                        "~/scratch/TCGA.SKCM.sumstats.tsv.gz"),
#              "name"=c("PROFILE", "TCGA"),
#              "outfile"="~/scratch/meta.test.tsv")

# Sanity-check lengths of --stats and --names
if(is.null(args$name)){
  args$name <- paste("stats", 1:length(args$stats), sep="")
}else{
  if(length(args$stats) != length(args$name)){
    stop("Different number of --stats and --name arguments provided. Exiting.")
  }
}
n.cohorts <- length(args$stats)

# Load each stats input
stats.list <- apply(data.frame(args$stats, args$name), 1,
                    function(rvals){
                      load.stats.single(rvals[1], suffix=rvals[2])
                    })
names(stats.list) <- args$name

# Merge stats across cohorts & conduct meta-analysis
stats <- merge.stats(stats.list)

# Write cleaned stats to --outfile
colnames(stats)[1] <- paste("#", colnames(stats)[1], sep="")
write.table(stats, args$outfile, col.names=T, row.names=F, sep="\t", quote=F)
