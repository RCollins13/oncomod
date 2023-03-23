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
require(bedr, quietly=TRUE)
require(GenomicRanges, quietly=TRUE)
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
parser$add_argument("--gtf", metavar=".gft", type="character",
                    help="GTF file for gene annotations")
parser$add_argument("--gw-sig", metavar="float", type="numeric", default=10e-8,
                    help="Specify strictest P-value threshold to annotate [default: 10e-8]")
parser$add_argument("--lenient-sig", metavar="float", type="numeric",
                    help=paste("If desired, specify a secondary, more lenient",
                               "P-value threshold to annotate in addition to",
                               "--gw-sig"))
parser$add_argument("--cancer-type", metavar="character",
                    help=paste("Color based on this cancer type"))
parser$add_argument("--highlight-gene", metavar="string", type="character",
                    help="Specify gene symbol to be highlighted, if optioned.")
parser$add_argument("--title", help="Custom plot title [default: --somatic-set-id]")
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="output .png file for plot")
args <- parser$parse_args()

# # DEV:
# args <- list("stats" = "~/scratch/pooled.CRAD.sumstats.tsv.gz",
#              # "somatic_set_id" = "COMUT|ENST00000256078_p.Gly12Val|KRAS_AMP",
#              # "somatic_set_id" = "COMUT|ENST00000256078_p.Gly12Asp|KRAS_AMP",
#              "somatic_set_id" = "ENST00000256078_p.Gly12Asp",
#              "gtf" = "~/scratch/gencode.v19.canonical.gtf.gz",
#              "gw_sig" = 0.0000006155892,
#              "lenient_sig" = 0.05/27202,
#              "cancer_type" = "CRAD",
#              "highlight_gene" = "KRAS",
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

# Load gene info, if optioned
if(!is.null(args$gtf)){
  genes <- load.gtf(args$gtf, paste(chrom, ":", xlims[1], "-", xlims[2], sep=""),
                    clean=TRUE)
  ylims[1] <- -0.15 * ylims[2]
  if(!is.null(args$highlight_gene)){
    if(any(genes$gene_name == args$highlight_gene)){
      h.gtf <- genes[genes$gene_name == args$highlight_gene]
      h.e.df <- data.frame(ranges(h.gtf[h.gtf$type == "exon"]))
      h.g.df <- data.frame(ranges(h.gtf[h.gtf$type == "gene"]))
      h.label <- args$highlight_gene
    }else{
      h.label <- NULL
    }
  }
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
if(!is.null(h.label)){
  rect(xleft=h.coords[1], xright=h.coords[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       bty="n", border=NA, col=adjustcolor(exon.color, alpha=0.2))
}

# Plot genes, if optioned
if(!is.null(args$gtf)){
  if(!is.null(h.label)){
    g.boring <- "gray30"
    e.boring <- "gray50"
  }else{
    g.boring <- gene.color
    e.boring <- exon.color
  }
  g.coords <- data.frame(ranges(genes[genes$type == "gene"]))
  segments(x0=g.coords[, 1], x1=g.coords[, 2], y0=0.5*ylims[1], y1=0.5*ylims[1],
           lwd=2, col=g.boring)
  e.coords <- data.frame(ranges(genes[genes$type == "exon"]))
  rect(xleft=e.coords[, 1], xright=e.coords[, 2],
       ybottom=0.75*ylims[1], ytop=0.25*ylims[1],
       col=e.boring, border=e.boring)
  if(!is.null(h.label)){
    segments(x0=h.g.df[, 1], x1=h.g.df[, 2], y0=0.5*ylims[1], y1=0.5*ylims[1],
             lwd=3, col=gene.color)
    rect(xleft=h.e.df[, 1], xright=h.e.df[, 2],
         ybottom=0.75*ylims[1], ytop=0.25*ylims[1],
         col=exon.color, border=exon.color, lwd=2)
  }
}

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
if(!is.null(args$gtf)){
  abline(h=0)
  axis(2, at=0.5*ylims[1], tick=F, line=-0.8, las=2, cex.axis=5/6,
       col.axis=g.boring, labels="Genes")
  if(!is.null(h.label)){
    axis(2, at=par("usr")[3], tick=F, line=-0.8, las=2, cex.axis=5/6,
         font=3, col.axis=exon.color, labels=h.label)
  }
}
mtext(3, font=2, text=plot.title)

# Close plot device
dev.off()

