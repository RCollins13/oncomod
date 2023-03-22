#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Generate a local Manhattan plot of germline association statistics


#########
# Setup #
#########
# Load necessary libraries and constants
require(RASMod, quietly=TRUE)
require(argparse, quietly=TRUE)
RASMod::load.constants("all")


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Plot germline association statistics",
                                           "for a locus of interest"))
parser$add_argument("--stats", metavar=".tsv", type="character",
                    help="association stats", required=TRUE)
parser$add_argument("--somatic-set-id", metavar="string", type="character",
                    help="somatic set identifier", required=TRUE)
parser$add_argument("--gw-sig", metavar="float", type="numeric", default=10e-8,
                    help="Specify strictest P-value threshold to annotate [default: 10e-8]")
parser$add_argument("--lenient-sig", metavar="float", type="numeric",
                    help=paste("If desired, specify a secondary, more lenient",
                               "P-value threshold to annotate in addition to",
                               "--gw-sig"))
parser$add_argument("--cancer-type", metavar="character",
                    help=paste("Color based on this cancer type"))
parser$add_argument("--title", help="Custom plot title [default: --somatic-set-id]")
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="output .png file for plot")
args <- parser$parse_args()

# # DEV:
# args <- list("stats" = "~/scratch/pooled.CRAD.sumstats.tsv.gz",
#              # "somatic_set_id" = "COMUT|ENST00000256078_p.Gly12Val|KRAS_AMP",
#              "somatic_set_id" = "COMUT|ENST00000256078_p.Gly12Asp|KRAS_AMP",
#              "gw_sig" = 0.0000006155892,
#              "lenient_sig" = 0.05/27202,
#              "cancer_type" = "CRAD",
#              "outfile" = "~/scratch/local_manhattan_test.png")

# Load stats input
stats <- load.assoc.stats.single(args$stats)

# Subset stats to somatic set of interest
stats <- stats[which(stats$somatic == args$somatic_set_id), ]

# Infer position of germline variants and only keep those with defined positions
stats$pos <- sapply(stats$germline, infer.germline.position.from.id)
stats <- stats[which(!is.na(stats$pos)), ]
stats$log.p <- -log10(stats$p)

# Infer plot parameters
xlims <- range(stats$pos, na.rm=T)
ylims <- c(0, 1.05 * max(stats$log.p, na.rm=T))
chrom <- names(sort(table(sapply(stats[1:20, "germline"],
                                        infer.germline.position.from.id,
                                        extract="chromosome")),
                           decreasing=T))[1]
plot.title <- if(is.null(args$title)){args$somatic_set_id}else{args$title}
if(is.null(args$cancer_type)){
  col.pal <- PANCAN.colors
}else{
  col.pal <- cancer.palettes[[args$cancer_type]]
}

# Add pointwise graphical parameters to stats
stats$pch <- sapply(stats$beta, function(b){if(b<0){25}else{24}})
stats$col <- col.pal[["light2"]]
stats$cex <- 0.3
if(!is.null(args$lenient_sig)){
  lenient.idxs <- which(stats$p <= args$lenient_sig)
  stats$col[lenient.idxs] <- col.pal[["main"]]
  stats$cex[lenient.idxs] <- 0.4
}
gw.idxs <- which(stats$p <= args$gw_sig)
stats$col[gw.idxs] <- col.pal[["dark2"]]
stats$cex[gw.idxs] <- 0.5

# Open plot device
png(args$outfile, height=3*350, width=5*350, res=350)

# Prep plot area
prep.plot.area(xlims, ylims, c(2.5, 2.5, 1, 0.1))

# Add significance thresholds
abline(h=-log10(c(args$gw_sig, if(!is.null(args$lenient_sig)){args$lenient_sig}else{NA})),
       col=col.pal[c("dark2", "main")], lty=c(5, 2))

# Add points
points(stats$pos, stats$log.p, pch=stats$pch, cex=stats$cex, bg=stats$col, col=NA, lwd=0)

# Add axes and title
clean.axis(1, at=axTicks(1), labels=paste(round(axTicks(1)/1000000, 1), "Mb", sep=""),
           infinite=T, title=paste("Position on", chrom))
clean.axis(2, at=unique(c(0, axTicks(2)[which(axTicks(2) > 0)])),
           title=bquote(-log[10](italic(P))))
mtext(3, font=2, text=plot.title)

# Close plot device
dev.off()

