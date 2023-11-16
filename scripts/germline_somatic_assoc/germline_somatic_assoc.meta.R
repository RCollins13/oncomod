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
require(metafor, quietly=TRUE)
require(argparse, quietly=TRUE)
require(OncoModR, quietly=TRUE)
OncoModR::load.constants("names")


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
parser$add_argument("--model", metavar="string", type="character", default="FE",
                    help="Specify meta-analysis model [Default: \"FE\"]")
parser$add_argument("--drop-frequencies", action="store_true", default=FALSE,
                    help="Do not report variant frequencies in --outfile")
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="output .tsv file for meta-analysis statistics")
args <- parser$parse_args()

# # DEV
# args <- list("stats"=c("~/scratch/PROFILE.PDAC.sumstats.tsv.gz",
#                        "~/scratch/TCGA.PDAC.sumstats.tsv.gz",
#                        "~/scratch/HMF.PDAC.sumstats.tsv.gz"),
#              "name"=c("PROFILE", "TCGA", "HMF"),
#              "model"="FE",
#              "drop_frequencies"=FALSE,
#              "outfile"="~/scratch/meta.test.tsv")

# # DEV: PAN-CANCER
# args <- list("stats" = c("~/scratch/PDAC.meta.sumstats.tsv.gz",
#                        "~/scratch/CRAD.meta.sumstats.tsv.gz",
#                        "~/scratch/LUAD.meta.sumstats.tsv.gz",
#                        "~/scratch/SKCM.meta.sumstats.tsv.gz"),
#              "name" = c("PDAC", "CRAD", "LUAD", "SKCM"),
#              "model" = "REML",
#              "drop_frequencies" = TRUE,
#              "outfile" = "~/scratch/meta.test.tsv")


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
                      load.assoc.stats.single(rvals[1], suffix=rvals[2])
                    })
names(stats.list) <- args$name

# Merge stats across cohorts
stats <- merge.assoc.stats(stats.list)

# Average variant frequencies, unless explicitly disabled
if(args$drop_frequencies){
  cols.to.drop <- grep("^samples\\.|^somatic_AC\\.|\\.germline_AC\\.|germline_AF\\.|^somatic_AF\\.",
                       colnames(stats))
  if(length(cols.to.drop) > 0){
    stats <- stats[, -cols.to.drop]
  }
}else{
  stats <- average.meta.freqs(stats)
}

# Meta-analyze effect sizes
stats <- germline.somatic.meta.analysis(stats, args$model)

# Add FDR Q-value
stats$fdr_q <- p.adjust(stats$p, method="fdr")

# Write cleaned stats to --outfile
colnames(stats)[1] <- paste("#", colnames(stats)[1], sep="")
write.table(stats, args$outfile, col.names=T, row.names=F, sep="\t", quote=F)
