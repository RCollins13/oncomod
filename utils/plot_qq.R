#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Generate a quantile-quantile (QQ) plot for a set of association statistics


#########
# Setup #
#########
# Load necessary libraries and constants
require(rCNV2, quietly=TRUE)
require(RASMod, quietly=TRUE)
require(argparse, quietly=TRUE)
RASMod::load.constants(c("names", "colors"))


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Generate a quantile-quantile (QQ)",
                                           "plot for a set of association statistics"))
parser$add_argument("--stats", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="output .png file for QQ plot")
parser$add_argument("--cancer", metavar="string",
                    help=paste("Color based on this cancer type [default: grey]"))
parser$add_argument("--cohort", metavar="string",
                     help=paste("Color based on this cohort [default: color by cancer]"))
parser$add_argument("--p-column", metavar="string", default="p",
                    help=paste("Column name for P-values [default: \"p\"]"))
parser$add_argument("--p-threshold", metavar="numeric", type="double",
                    help=paste("P-value threshold for significance [default: Bonferroni]"))
parser$add_argument("--pt-cex", type="double", metavar="float", default=0.35,
                    help=paste("Size scalar for points [default: 0.35]"))
parser$add_argument("--smallest-p", type="double", metavar="float", default=1e-20,
                    help=paste("P-values below this threshold will be rounded [default: 10E-20]"))
parser$add_argument("--title", metavar="string", type="character", help="Custom plot title")
args <- parser$parse_args()

# # DEV
# args <- list("stats" = "~/scratch/TCGA.LUAD.sumstats.tsv.gz",
#           "outfile" = "~/scratch/TCGA.LUAD.qq.png",
#           "cancer" = "LUAD",
#           "cohort" = "TCGA",
#           "p_column" = "p",
#           "p_threshold" = 0.0000003218497,
#           "pt_cex" = 0.35,
#           "smallest_p" = 10E-20)

# Load data
df <- read.table(args$stats, header=T, sep="\t", check.names=F, comment.char="")
if(args$p_column != "p"){
  colnames(df)[which(colnames(df) == args$p_column)] <- "p"
}

# Determine coloring
if(!is.null(args$cancer)){
  title <- cancer.names.short[args$cancer]
  if(!is.null(args$cohort)){
    if(args$cohort %in% names(cohort.names.short)){
      title <- paste(title, " (", cohort.names.short[args$cohort], ")", sep="")
      color <- cancer.palettes[[args$cancer]][cohort.color.prefixes[args$cohort]]
    }else{
      title <- paste(title, " (", args$cohort, ")", sep="")
      color <- cancer.colors[args$cancer]
    }
  }else{
    color <- cancer.colors[args$cancer]
  }
}else{
  color <- "gray15"
  title <- NULL
}
if(!is.null(args$title)){
  title <- args$title
}

# Plot QQ
png(args$outfile, height=3*300, width=3*300, res=300)
rCNV2::plot.qq(df, smallest.p=args$smallest_p, cutoff=args$p_threshold,
               pt.color=color, pt.cex=args$pt_cex, echo.lambdas=T)
mtext(title, 3, line=-0.5, cex=0.85)
dev.off()
