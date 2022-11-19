#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Generate a series of descriptive plots for cohort- and cancer-level phenotype metadata


#########
# Setup #
#########
# Load necessary libraries and constants
require(RASMod, quietly=TRUE)
require(argparse, quietly=TRUE)
RASMod::load.constants("all")


##################
# Data functions #
##################


######################
# Plotting functions #
######################
# PC1 x PC2 with stacked bar of populations for single cohort
pc.scatter.w.bar <- function(df, colors, pt.cex=0.25, bar.wex=0.15,
                             parmar=c(2.5, 3, 0.25, 3.5)){
  # Get plot values
  x <- as.numeric(df$PC1)
  pc.xlims <- range(x, na.rm=T)
  bar.spacer <- 0.5 * bar.wex * diff(pc.xlims)
  xlims <- c(pc.xlims[1], pc.xlims[2] + (3 * bar.spacer))
  y <- as.numeric(df$PC2)
  ylims <- range(y, na.rm=T)
  pops <- df$POPULATION
  pw.colors <- colors[pops]

  # Get mean Y position and proportion per population to determine bar/legend order
  pop.y <- sort(sapply(unique(pops), function(pop){
    mean(y[which(pops == pop)], na.rm=T)
  }))
  pop.order <- names(sort(pop.y))
  pop.counts <- table(pops)[pop.order]
  pop.pct <- pop.counts / sum(pop.counts)
  pop.pct.cumul <- cumsum(pop.pct)
  pop.pct.scaled <- (pop.pct.cumul * diff(ylims)) + ylims[1]

  # Prep plot area
  RASMod::prep.plot.area(xlims, ylims, parmar, xaxs="i", yaxs="i")

  # Add points
  points(x, y, cex=pt.cex, col=pw.colors, pch=19, xpd=T)

  # Add axes
  axis(1, at=c(-10e10, pc.xlims[2]), tck=0, labels=NA)
  RASMod::clean.axis(1, at=axTicks(1)[which(axTicks(1) < pc.xlims[2])], title="PC1")
  RASMod::clean.axis(2, infinite=TRUE, title="PC2")

  # Add legend
  legend.y.at <- (c(ylims[1], pop.pct.scaled[-length(pop.order)]) + pop.pct.scaled) / 2
  yaxis.legend(pop.order, x=xlims[2], legend.y.at, sep.wex=0.05 * diff(xlims),
               min.label.spacing=0.075 * diff(ylims),
               lower.limit=ylims[1] + (0.025 * diff(ylims)),
               upper.limit=ylims[1] + (0.975 * diff(ylims)),
               colors=colors[pop.order])

  # Add bars corresponding to population representation
  rect(xleft=pc.xlims[2] + bar.spacer, xright=xlims[2],
       ybottom=c(ylims[1], pop.pct.scaled[-length(pop.order)]),
       ytop=pop.pct.scaled, border=NA, bty="n", col=colors[pop.order])
  rect(xleft=pc.xlims[2] + bar.spacer, xright=xlims[2],
       ybottom=ylims[1], ytop=ylims[2], col=NA, xpd=TRUE)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Summarize cohort and patient phenotype metadata")
parser$add_argument('--metadata', metavar='.tsv', type="character", nargs='+',
                    help='cohort metadata .tsv', action='append', required=TRUE)
parser$add_argument('--cohort-name', metavar='string', type="character", nargs='*',
                    action='append', help='cohort names (order must batch --metadata)')
parser$add_argument('--out-prefix', metavar='path', type="character", nargs=1,
                    help='file prefix for all plots')
args <- parser$parse_args()
meta.in <- args$metadata[, 1]
if(is.null(args$cohort_name)){
  names(meta.in) <- paste("Cohort", 1:length(meta.in))
}else{
  names(meta.in) <- args$cohort_name
}
out.prefix <- args$out_prefix

# # Dev
# meta.in <- c("TCGA" = "~/scratch/TCGA.ALL.sample_metadata.tsv.gz",
#              "DFCI" = "~/scratch/PROFILE.ALL.sample_metadata.tsv.gz")
# out.prefix <- "~/scratch/rasmod.test"

# Load metadata
meta <- lapply(meta.in, load.patient.metadata)


#####################
# Cohort-wide plots #
#####################
# Plot distribution of patients by cancer type per cohort
pdf(paste(out.prefix, "cancer_type_by_cohort.bars.pdf", sep="."),
    height=3, width=3)
scaled.bars(lapply(meta, function(df){df$CANCER_TYPE}),
            colors=cancer.colors, title="Cancer")
dev.off()


# Plot distribution of patients by sex per cohort
pdf(paste(out.prefix, "sex_by_cohort.bars.pdf", sep="."),
    height=3, width=3)
scaled.bars(lapply(meta, function(df){df$SEX}),
            colors=sex.colors, title="Sex")
dev.off()


#########################
# Cohort-specific plots #
#########################
for(cohort in names(meta)){
  # Scatterplot of PC1xPC2 colored by ancestry
  pdf(paste(out.prefix, cohort, "ancestry_pcs.pdf", sep="."),
      height=2.5, width=3.25)
  pc.scatter.w.bar(meta[[cohort]], pop.colors)
  dev.off()
}


#########################
# Cancer-specific plots #
#########################
cancer.types <- unique(unlist(sapply(meta, function(df){df$CANCER_TYPE})))
sapply(cancer.types, function(cancer){
  # Subset patient metadata to only those patients from cancer type
  meta.sub <- lapply(meta, function(df){df[which(df$CANCER_TYPE == cancer), ]})

  # Make subdirectory for plots from this cancer type, if necessary
  sub.prefix <- paste(dirname(out.prefix), cancer, basename(out.prefix), sep="/")
  if(!dir.exists(dirname(sub.prefix))){
    dir.create(dirname(sub.prefix))
  }

  # Plot distribution of patients by ancestry per cohort
  pdf(paste(sub.prefix, cancer, "ancestry_by_cancer_type.bars.pdf", sep="."),
      height=3, width=3)
  scaled.bars(lapply(meta.sub, function(df){df$POPULATION}),
              colors=pop.colors, title="Population")
  dev.off()

  # Plot distribution of patients by sex per cohort
  pdf(paste(sub.prefix, cancer, "sex_by_cohort.bars.pdf", sep="."),
      height=3, width=3)
  scaled.bars(lapply(meta.sub, function(df){df$SEX}),
              colors=sex.colors, title="Sex")
  dev.off()

  # Plot distribution of patients by stage per cohort
  pdf(paste(sub.prefix, cancer, "stage_by_cohort.bars.pdf", sep="."),
      height=3, width=3)
  scaled.bars(lapply(meta.sub, function(df){df$APPROX_STAGE}),
              colors=stage.colors, legend.names=stage.names, title="Stage")
  dev.off()
})

