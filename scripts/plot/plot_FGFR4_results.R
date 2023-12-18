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
require(beeswarm, quietly=TRUE)
load.constants("all")

# Two simple positional arguments
args <- commandArgs(trailingOnly=TRUE)
tsv.in <- as.character(args[1])
dir.out <- as.character(args[2])

# # DEV:
# tsv.in <- "~/scratch/pooled.CRAD.FGFR4_allelic_series.sumstats.tsv"
# dir.out <- "~/scratch/"

# Read FGFR4:KRAS association data
data <- read.table(tsv.in, header=T, sep="\t", comment.char="", check.names=F)
colnames(data)[1] <- gsub("^#", "", colnames(data)[1])

# Infer missing data
data$AF <- (data$yes_somatic.germline_AC + data$no_somatic.germline_AC) / (2 * data$samples)
data$somatic_freq <- data$somatic_AC / data$samples
data$csq <- "missense"
data$csq[grep("=$", data$germline)] <- "synonymous"
data$csq[grep("Ter$|Ter[0-9]+$", data$germline)] <- "lof"
data$aa <- sapply(data$germline, function(gstr){
  as.numeric(sub("=", "", unlist(strsplit(unlist(strsplit(gstr, split="_p.", fixed=T))[2], split="[A-z]+"))[2]))
})

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

# Map FGFR4 variants onto protein domains by frequency
protein.regions <- list("Extracellular" = c(22, 369),
                        "Transmem." = c(370, 390),
                        "Cytoplasmic" = c(391, 802))
disulfide.bonds <- list("DSB 1" = c(57, 101),
                        "DSB 2" = c(172, 224),
                        "DSB 3" = c(271, 333))
protein.domains <- list("Ig-C2 1" = c(22, 118),
                        "Disordered" = c(119, 148),
                        "Ig-C2 2" = c(152, 240),
                        "Ig-C2 3" = c(249, 349),
                        "Kinase" = c(467, 755))
phosphosites <- c(390, 419, 573, 642, 643, 754)
pdf(paste(dir.out, "FGFR4.annotated_variant_localization.pdf", sep="/"),
    height=4.2, width=6)
par(mar=c(0.25, 5, 2.5, 0.25), bty="n")
plot(NA, xlim=c(0, 802), ylim=c(7, 0),
     xaxt="n", yaxt="n", xlab="", ylab="")
rect(xleft=c(unlist(disulfide.bonds)[c(3, 5)],
             protein.regions[[2]][1]),
     xright=c(unlist(disulfide.bonds)[c(4, 6)],
             protein.regions[[2]][2]),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     border=NA, bty="n", col=adjustcolor("yellow", 0.3))
abline(h=3:7)
sapply(1:length(protein.regions), function(i){
  rect(xleft=protein.regions[[i]][1], xright=protein.regions[[i]][2],
       ybottom=0.1, ytop=0.9, col="azure3")
})
sapply(1:length(protein.regions), function(i){
  text(x=mean(protein.regions[[i]]), y=0.5, labels=names(protein.regions)[i],
       cex=0.55, srt=45)
})
sapply(1:length(disulfide.bonds), function(i){
  rect(xleft=disulfide.bonds[[i]][1], xright=disulfide.bonds[[i]][2],
       ybottom=1.1, ytop=1.9, col="cadetblue")
})
sapply(1:length(disulfide.bonds), function(i){
  text(x=mean(disulfide.bonds[[i]]), y=1.5, labels="Disulfide\nBond", cex=0.55, srt=45)
})
sapply(1:length(protein.domains), function(i){
  rect(xleft=protein.domains[[i]][1], xright=protein.domains[[i]][2],
       ybottom=2.1, ytop=2.9, col="lightcyan3")
})
sapply(1:length(protein.domains), function(i){
  text(x=mean(protein.domains[[i]]), y=2.5, labels=names(protein.domains)[i],
       cex=0.55, srt=45)
})
common.vars <- t1.data[which(t1.data$AF >= 0.01), ]
cv <- beeswarm(common.vars$aa, at=3.5, add=T, horiz=T, pch=19,
               pwcol=csq.colors[common.vars$csq], corral="wrap", corralWidth=0.8)
common.labels <- label.vars[sapply(cv$y.orig, function(y){t1.data$germline[which(t1.data$aa == y & t1.data$AF >= 0.01)]})]
common.nonna <- which(!is.na(common.labels))
common.labels <- common.labels[common.nonna]
common.colors <- csq.colors[t1.data$csq[unlist(sapply(names(common.labels), function(sid){which(t1.data$germline == sid)}))]]
text(x=cv$y.orig[common.nonna], y=3.5, pos=3, labels=common.labels, cex=4/6, col=common.colors)
rare.case.vars <- t1.data[which(t1.data$AF < 0.01
                                & t1.data$no_somatic.germline_AC == 0), ]
beeswarm(rare.case.vars$aa, at=4.5, add=T, horiz=T, pch=19,
         pwcol=csq.colors[rare.case.vars$csq], corral="wrap", corralWidth=0.8)
rare.wt.vars <- t1.data[which(t1.data$AF < 0.01
                              & t1.data$yes_somatic.germline_AC == 0), ]
beeswarm(rare.wt.vars$aa, at=5.5, add=T, horiz=T, pch=19,
         pwcol=csq.colors[rare.wt.vars$csq], corral="wrap", corralWidth=0.8)
rare.both.vars <- t1.data[which(t1.data$AF < 0.01
                                & t1.data$yes_somatic.germline_AC > 0
                                & t1.data$no_somatic.germline_AC > 0), ]
beeswarm(rare.both.vars$aa, at=6.5, add=T, horiz=T, pch=19,
         pwcol=csq.colors[rare.both.vars$csq], corral="wrap", corralWidth=0.8)
clean.axis(3, at=seq(0, 1000, 100), title=bquote(italic("FGFR4") ~ "Residue"))
axis(2, at=(1:7)-0.5, las=2, cex.axis=5/6, tick=F, line=-1,
     labels=c("Localization", "Bonds", "Domains", "Common",
              "Rare vars in\nKRAS mutant", "Rare vars in\nKRAS WT",
              "Rare vars in\nmutant & WT"))
dev.off()

# Compute relative frequency of rare FGFR4 germline hits by domain
max.n.case <- max(t1.data$somatic_AC)
max.n.control <- max(t1.data$samples) - max.n.case
get.region.stats <- function(regions){
  do.call("rbind", lapply(1:length(regions), function(i){
    hits <- which(t1.data$aa >= regions[[i]][1]
                  & t1.data$aa <= regions[[i]][2]
                  & t1.data$csq != "synonymous"
                  & t1.data$AF < 0.01)
    case.hits <- length(intersect(which(t1.data$yes_somatic.germline_AC > 0), hits))
    control.hits <- length(intersect(which(t1.data$no_somatic.germline_AC > 0), hits))
    fisher.res <- unlist(fisher.test(matrix(c(max.n.control - control.hits,
                                              max.n.case - case.hits,
                                              control.hits, case.hits),
                                            byrow=T, nrow=2))[c("estimate", "p.value")])
    return(c(names(regions)[i], case.hits, control.hits, fisher.res))
  }))
}
region.burden <- as.data.frame(do.call("rbind", lapply(list(protein.regions, protein.domains, disulfide.bonds), get.region.stats)))
colnames(region.burden) <- c("region", "n.mut_only.vars", "n.WT_only.vars", "OR", "P")
write.table(region.burden[order(region.burden[, 4], decreasing=TRUE), ],
            paste(dir.out, "FGFR4.regional_burden.quick_and_dirty.tsv", sep="/"),
            col.names=T, row.names=F, sep="\t", quote=F)
