#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Generate summary plots of FGFR4 : KRAS single-variant associations


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(OncoModR, quietly=TRUE)
load.constants("all")

# Two simple positional arguments
args <- commandArgs(trailingOnly=TRUE)
tsv.in <- as.character(args[1])
dir.out <- as.character(args[2])

# Read FGFR4:KRAS association data
data <- read.table(tsv.in, header=T, sep="\t", comment.char="", check.names=F)
colnames(data)[1] <- gsub("^#", "", colnames(data)[1])

# Infer missing data
data$AF <- (data$yes_somatic.germline_AC + data$no_somatic.germline_AC) / (2 * data$samples)
data$somatic_freq <- data$somatic_AC / data$samples
data$csq <- "missense"
data$csq[grep("=$", data$germline)] <- "synonymous"
data$csq[grep("Ter$|Ter[0-9]+$", data$germline)] <- "lof"

# Custom labels
label.vars <- c("ENST00000292408_p.Gly388Arg" = "G388R",
  "ENST00000292408_p.Pro136Leu" = "P136L",
  "ENST00000292408_p.Arg54=" = "R54R",
  "ENST00000292408_p.Ser137=" = "S137S",
  "ENST00000292408_p.Asp575Asn" = "D575N",
  "ENST00000292408_p.Val10Ile" = "V10I",
  "ENST00000292408_p.Arg234=" = "R234R")

# Stacked barplot of coding variants by AF & consequence
af.breaks <- seq(0,1, 0.05)
bdat <- do.call("cbind", lapply(names(csq.colors), function(csq){
  hist(data$AF[intersect(which(data$csq == csq), grep("tier1", data$somatic))],
       breaks=af.breaks, plot=F)$counts
}))
pdf(paste(dir.out, "FGFR4_alleles.AF_hist.pdf", sep="/"),
    height=3.25, width=3.25)
par(mar=c(2.5, 3, 2, 0.5), bty="n")
plot(NA, xlim=c(0, 1), ylim=range(apply(bdat, 1, sum)),
     xlab="", xaxt="n", ylab="", yaxt="n", xaxs="i", yaxs="i")
rect(xleft=af.breaks[-length(af.breaks)], xright=af.breaks[-1],
     ybottom=0, ytop=bdat[, 1], col=csq.colors["synonymous"])
rect(xleft=af.breaks[-length(af.breaks)], xright=af.breaks[-1],
     ybottom=bdat[, 1], ytop=apply(bdat[, 1:2], 1, sum),
     col=csq.colors["missense"])
rect(xleft=af.breaks[-length(af.breaks)], xright=af.breaks[-1],
     ybottom=apply(bdat[, 1:2], 1, sum), ytop=apply(bdat[, 1:3], 1, sum),
     col=csq.colors["lof"])
clean.axis(1, title="Allele Frequency", infinite=TRUE)
clean.axis(2, title="Count", infinite=TRUE)
mtext(3, text=bquote(italic("FGFR4") ~ "Variants by Frequency"), line=0.5)
legend("topright", cex=5/6, fill=csq.colors, bty="n",
       legend=c("Synonymous", "Missense", "Loss-of-Function"))
dev.off()

# Volcano plot of KRAS tier 1
t1.data <- data[grep("tier1", data$somatic), ]
pdf(paste(dir.out, "KRAS_tier_1.FGFR4_alleles.volcano.pdf", sep="/"),
    height=3.25, width=3.25)
par(mar=c(2.5, 3, 2, 0.5), bty="n")
plot(log2(exp(t1.data$beta)), -log10(t1.data$p), pch=19, col=csq.colors[t1.data$csq],
     xlab="", xaxt="n", ylab="", yaxt="n",
     panel.first=c(abline(h=-log10(0.05), v=0, lty=5, col="gray60")))
sapply(names(label.vars), function(var){
  text(x=log2(exp(t1.data$beta[which(t1.data$germline == var)])),
  y=-log10(t1.data$p[which(t1.data$germline == var)]),
  labels=label.vars[var], pos=2, cex=5/6,
  col=csq.colors[t1.data$csq[which(t1.data$germline == var)]])
})
clean.axis(1, title=bquote(log[2]("Odds Ratio")), infinite=TRUE)
clean.axis(2, title=bquote(-log[10](italic(P))), infinite=TRUE)
mtext(3, font=2, text=bquote(italic(FGFR4) ~ "germline alleles vs."), line=1)
mtext(3, font=2, text=bquote(italic(KRAS) ~ "tier 1 somatic mutations"), line=0)
legend("topright", cex=5/6, pch=19, col=csq.colors, legend=csq.names.short)
dev.off()

# Scatterplot of tier 1 association beta vs. allele frequency
pdf(paste(dir.out, "KRAS_tier_1.FGFR4_alleles.beta_vs_AF.pdf", sep="/"),
    height=3.25, width=3.25)
par(mar=c(2.5, 3, 2, 0.5), bty="n")
plot(log2(exp(t1.data$beta)), t1.data$AF, pch=19, col=csq.colors[t1.data$csq],
     xlab="", xaxt="n", ylab="", yaxt="n",
     panel.first=c(abline(v=0, lty=5, col="gray60")))
sapply(names(label.vars), function(var){
  text(x=log2(exp(t1.data$beta[which(t1.data$germline == var)])),
       y=t1.data$AF[which(t1.data$germline == var)],
       labels=label.vars[var], pos=2, cex=5/6,
       col=csq.colors[t1.data$csq[which(t1.data$germline == var)]])
})
clean.axis(1, title=bquote(log[2]("Odds Ratio")), infinite=TRUE)
clean.axis(2, title="Allele Frequency", infinite=TRUE)
mtext(3, font=2, text=bquote(italic(FGFR4) ~ "germline alleles vs."), line=1)
mtext(3, font=2, text=bquote(italic(KRAS) ~ "tier 1 somatic mutations"), line=0)
legend("topright", cex=5/6, pch=19, col=csq.colors, legend=csq.names.short)
dev.off()

