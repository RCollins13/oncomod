#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Generate full-slide summary diagram of KRAS somatic mutations in dataset


#########
# Setup #
#########
# Load necessary libraries and constants
require(OncoModR, quietly=TRUE)
require(argparse, quietly=TRUE)
require(meta, quietly=TRUE)
require(viridis, quietly=TRUE)
require(shape, quietly=TRUE)
OncoModR::load.constants("all")


##################
# Data functions #
##################
# Load & clean somatic frequency data from .tsv
load.somatic.freqs <- function(tsv.in){
  # Read data into memory
  data <- read.table(tsv.in, sep="\t", comment.char="", check.names=F, header=T)
  rownames(data) <- data$csq
  data <- data[which(data$gene == "KRAS"), ]

  # Drop unnecessary columns
  drop.cols <- c("position", "gene",
                 colnames(data)[grep("_set_id$|vids$", colnames(data))])
  data <- data[, -which(colnames(data) %in% drop.cols)]

  # Add new row for combined frequency of any missense
  cohorts <- unique(unlist(sapply(strsplit(colnames(data)[grep("_AF$", colnames(data))],
                                           split="_"), function(l){l[1]})))
  cancers <- unique(unlist(sapply(strsplit(colnames(data)[grep("_AF$", colnames(data))],
                                           split="_"), function(l){l[2]})))
  mis.idxs <- which(data$alt %in% toupper(letters))
  all.mis <- do.call("cbind", lapply(cohorts, function(cohort){
    cohort.res <- do.call("cbind", sapply(cancers, function(cancer){
      max.an <- max(data[mis.idxs, paste(cohort, cancer, "AN", sep="_")], na.rm=T)
      sum.ac <- sum(data[mis.idxs, paste(cohort, cancer, "AC", sep="_")], na.rm=T)
      af <- sum.ac / max.an
      data.frame(max.an, sum.ac, af)
    }))
    colnames(cohort.res) <- as.vector(sapply(cancers, function(cancer){
      paste(cohort, cancer, c("AN", "AC", "AF"), sep="_")
    }))
    return(cohort.res)
  }))
  other.mis.info <- data.frame("chrom" = NA, "csq" = "any_missense", "codon" = NA,
                               "ref" = NA, alt = NA)
  all.mis <- cbind(other.mis.info, all.mis)
  data <- rbind(data, all.mis)
  rownames(data)[nrow(data)] <- "any_mis"

  # Compute inverse variance-weighted mean freq for all mutations per cancer type
  for(cancer in cancers){
    mfreqs <- as.data.frame(do.call("rbind", lapply(1:nrow(data), function(i){
      n.all <- as.numeric(data[i, grep(paste("", cancer, "AN", sep="_"), colnames(data))])
      n.hit <- as.numeric(data[i, grep(paste("", cancer, "AC", sep="_"), colnames(data))])
      meta.res <- metaprop(n.hit, n.all, method="Inverse")
      meta:::backtransf(unlist(meta.res[c("TE.random", "lower.random", "upper.random")]),
                        sm=meta.res$sm)
    })))
    colnames(mfreqs) <- paste("meta_", cancer, "_AF", c("", "_lower", "_upper"), sep="")
    data <- cbind(data, mfreqs)
  }

  # Return sorted df
  data[, sort(colnames(data))]
}

# Get count of unique alleles observed per codon
get.codon.mut.counts <- function(data, n.codons=189, max.hits=5){
  counts <- sapply(1:n.codons, function(k){length(which(data$codon == k))})
  counts[which(counts > max.hits)] <- max.hits
  return(counts)
}

# Compute IVW meta-analysis of fraction of patients with a mutation per residue
get.codon.mut.freq <- function(data, cancer, n.codons=189, max.freq=0.05,
                               scale.factor=10){
  cohorts <- unique(cohort.names.short)
  freq <- sapply(1:n.codons, function(k){
    k.idx <- which(data$codon == k)
    if(length(k.idx) == 0){return(0)}
    an <- apply(data[k.idx, grep(paste(cancer, "AN", sep="_"), colnames(data))],
                2, max, na.rm=T)
    ac <- apply(data[k.idx, grep(paste(cancer, "AC", sep="_"), colnames(data))],
                2, sum, na.rm=T)
    meta.res <- metaprop(ac, an, method="Inverse")
    meta:::backtransf(meta.res$TE.random, sm=meta.res$sm)
  })
  freq[which(freq > max.freq)] <- max.freq
  ceiling(scale.factor * 100 * freq)
}


######################
# Plotting functions #
######################
# Simple idiogram
kras.idio <- function(chrom.len=133851895, text.cex=5/6,
                      parmar=c(2, 3, 2, 1)){
  prep.plot.area(xlim=c(0, chrom.len), ylim=c(0, 1), parmar=parmar)
  segments(x0=0, x1=chrom.len, y0=0.5, y1=0.5)
  rect(xleft=25346761, xright=25415273, ybottom=0.3, ytop=0.7,
       col="red", border="red")
  axis(2, at=0.5, tick=F, line=-0.5, las=2, labels="chr12", cex.axis=text.cex)
}

# Simple gene structure plot
kras.gene <- function(view.start=25404589, view.stop=25357497,
                      exon.height=0.75, utr.height=0.3, gene.lwd=1.5,
                      text.cex=5/6, parmar=c(0.5, 2, 1, 1)){

  # Manually curate KRAS isoform coordinates from Gencode
  kras.4a.coords <- list("utr.starts" = c(25362365, 25368371, 25398319, 25403685),
                         "utr.stops" = c(25362845, 25368377, 25398329, 25403737),
                         "exon.starts" = c(25368378, 25378548, 25380168, 25398208),
                         "exon.stops" = c(25368494, 25378707, 25380346, 25398318))
  kras.4b.coords <- list("utr.starts" = c(25357723, 25398319, 25403685),
                         "utr.stops" = c(25362731, 25398329, 25403865),
                         "exon.starts" = c(25362732, 25378548, 25380168, 25398208),
                         "exon.stops" = c(25362845, 25378707, 25380346, 25398318))
  kras.coords <- list("4a" = kras.4a.coords, "4b" = kras.4b.coords)

  # Prepare plot layout
  prep.plot.area(xlim=c(view.start, view.stop), ylim=c(2.75, -0.75), parmar=parmar)
  clean.axis(3, labels=paste(axTicks(3) / 1000000, "Mb", sep=""), infinite=TRUE,
             tck=-0.05, col.axis="gray65", cex.axis=text.cex)
  text(x=max(unlist(unlist(kras.coords))), y=1, labels="KRAS", font=3,
       pos=2, xpd=T, cex=text.cex)

  # Add 4a & 4b coordinates
  segments(x0=sapply(kras.coords, function(l){min(unlist(l), na.rm=T)}),
           x1=sapply(kras.coords, function(l){max(unlist(l), na.rm=T)}),
           y0=c(0.5, 1.5), y1=c(0.5, 1.5), col=gene.color, lwd=gene.lwd)
  tss <- max(unlist(kras.coords[[1]]))
  segments(x0=tss, x1=tss, y=1-0.5, y1=1-0.65-(exon.height/2), col="gray65")
  Arrows(x0=tss, x1=tss+0.02*diff(par("usr")[1:2]),
         y0=1-0.65-(exon.height/2), y1=1-0.65-(exon.height/2),
         col="gray65", arr.type="triangle", arr.length=0.15, arr.width=0.075)
  sapply(1:2, function(i){
    rect(xleft=kras.coords[[i]]$utr.starts, xright=kras.coords[[i]]$utr.stops,
         ybottom=i-0.5-(utr.height/2), ytop=i-0.5+(utr.height/2),
         col=gene.color, border=gene.color)
    rect(xleft=kras.coords[[i]]$exon.starts, xright=kras.coords[[i]]$exon.stops,
         ybottom=i-0.5-(exon.height/2), ytop=i-0.5+(exon.height/2),
         col=gene.color, border=gene.color)
  })

  # Add labels for 4A/4B exons
  text(x=(kras.4a.coords$exon.starts + kras.4a.coords$exon.stops)/2,
       y=0.5-(utr.height/2), pos=3, labels=rev(c(1, 2, 3, "4A")),
       cex=text.cex, col=gene.color, xpd=T)
  text(x=((kras.4b.coords$exon.starts + kras.4b.coords$exon.stops)/2)[1],
       y=1.5+(utr.height/2), pos=1, labels="4B", cex=text.cex,
       col=gene.color, xpd=T)
}

# Draw 1d protein model of KRAS
kras.protein <- function(n.codons=189, label.domains=TRUE,
                         text.cex=5/6, parmar=c(0.5, 2, 0.5, 1)){
  # Features curated from literature & UniProt
  features <- data.frame("start"=c(10, 10, 29, 32, 30, 59, 60, 116, 166),
                         "stop"=c(18, 18, 35, 40, 38, 60, 76, 119, 189),
                         "class"=c("GTP", "P-Loop", "GTP", "Effector", "Switch I",
                                   "GTP", "Switch II", "GTP", "Hypervariable"),
                         "height"=c(1, 0.5, 1, 1, 0.5, 1, 0.5, 1, 1))
  exon.breaks <- c(37, 97, 150)

  # Set plot area
  prep.plot.area(xlims=c(0, n.codons), ylims=c(0, 1), parmar=parmar)

  # Draw rectangle
  rect(xleft=0, xright=n.codons, ybottom=0, ytop=1, xpd=T, border=NA, col="#EBE4DD")

  # Add domains
  rect(xleft=features$start, xright=features$stop,
       ybottom=0.5-(features$height/2), ytop=0.5+(features$height/2),
       border=protein.feat.colors[features$class],
       col=protein.feat.colors[features$class])

  # Add breaks for exons
  abline(v=exon.breaks, col="gray65")
  sapply(exon.breaks, function(x){
    axis(3, at=x, labels=NA, col.ticks="gray65", col="gray65", tck=-0.015)
  })

  # Add labels
  axis(2, at=0.5, tick=F, line=-0.7, labels="Protein\nDomains", cex.axis=text.cex, las=2)
  if(label.domains){
    text(x=apply(features[grep("Switch", features$class), c("start", "stop")], 1, mean),
         y=0.5, cex=text.cex, labels=c("Sw. I", "Switch II"))
    text(x=((features$start + features$stop)/2)[nrow(features)],
         y=0.5, cex=text.cex, col="white", labels="Hypervariable")
    text(x=((features$start + features$stop)/2)[which(features$class == "P-Loop")],
         y=0.5, cex=text.cex, labels="P-Loop")
  }
  axis(3, at=0, labels="N", line=-0.9, tick=F, col.axis="gray65", cex.axis=text.cex)
  axis(3, at=n.codons, labels="C", line=-0.9, tick=F, col.axis="gray65", cex.axis=text.cex)
}

# Simple legend for protein domains
protein.legend <- function(cex=5/6, parmar=c(0.1, 1.5, 0.1, 1.5)){
  prep.plot.area(xlims=c(0, 1), ylims=c(0, 1), parmar=parmar)
  legend("center",
         fill=protein.feat.colors[c("GTP", "Effector")],
         border=protein.feat.colors[c("GTP", "Effector")],
         legend=c("GTP Binding Site", "Effector Binding Site"),
         bg=NA, box.col=NA, box.lwd=0, cex=cex, xpd=T)
}

# Function for codon-by-codon heatmap of KRAS residues
codon.heat <- function(counts, pal, n.codons=189, text.cex=5/6,
                       left.label=NULL, number.residues=TRUE,
                       parmar=c(2, 2, 0.5, 0.5)){
  # Plot heatmap
  prep.plot.area(xlim=c(0, n.codons), ylim=c(0, 1), parmar=parmar)
  rect(xleft=(1:n.codons)-1, xright=1:n.codons, ybottom=0, ytop=1,
       col=pal[counts+1], border="gray65", xpd=T, lwd=0.5)
  rect(xleft=0, xright=n.codons, ybottom=0, ytop=1)
  axis(2, at=0.5, tick=F, las=2, labels=left.label, line=-0.9, cex.axis=text.cex)
  if(number.residues){
    axis(3, at=c(seq(0, n.codons, 40), n.codons), tick=F, line=-1,
         cex.axis=text.cex, col.axis="gray65")
  }
}

# Simple legend for codon allele heatmap
codon.mut.heat.legend <- function(max.hits=5, text.cex=5/6,
                                  parmar=c(0.1, 0.1, 0.1, 0.1)){
  prep.plot.area(xlims=c(0, max.hits+1), ylim=c(0, 2), parmar=parmar)
  pal <- colorRampPalette(stage.colors)(max.hits+1)
  rect(xleft=0:max.hits, xright=1:(max.hits+1), ybottom=0, ytop=1,
       col=pal[1:(max.hits+1)], xpd=T, border="gray65")
  text(x=(1:(max.hits))-0.5, y=0.5, labels=(1:(max.hits))-1, cex=text.cex)
  text(x=max.hits+0.5, y=0.5, labels=bquote("" >= .(max.hits)), cex=text.cex)
  text(x=mean(par("usr")[1:2]), y=1.5, labels="Unique Alleles at Site", xpd=T,
       cex=text.cex)
}

# Simple legend for codon mut freq heatmap
codon.mut.freq.heat.legend <- function(max.freq=0.05, text.cex=5/6,
                                       parmar=c(1.5, 2, 1.5, 2)){
  prep.plot.area(xlims=c(0, 100), ylim=c(0, 2), parmar=parmar)
  pal <- viridis(100)
  rect(xleft=0:99, xright=1:100, ybottom=0, ytop=1, col=pal, xpd=T, border=pal)
  rect(xleft=0, xright=100, ybottom=0, ytop=1, xpd=T, border="gray65")
  text(x=mean(par("usr")[1:2]), y=1.5, labels="Tumors with Mutation", xpd=T, cex=text.cex)
  axis(2, at=0.5, tick=F, line=-0.9, cex.axis=text.cex, las=2, labels="0%", xpd=T)
  axis(4, at=0.5, tick=F, line=-0.9, cex.axis=text.cex, las=2,
       labels=bquote("" >= .(paste(100*max.freq, "%", sep=""))), xpd=T)
}

# Helper function to plot a single panel of variant frequency meta-analyses
plot.freq.meta.single <- function(data, csq, cancer, title=NULL, text.cex=5/6,
                                  title.cex=0.65, diamond.ratio=1/3,
                                  proportional.bars=TRUE, x.space=1/10,
                                  ylims=NULL, y.axis=TRUE, label.y.axis=TRUE,
                                  title.y.axis=TRUE, title.line=0.75, ax.tck=-0.05,
                                  parmar=c(0.5, 3, 1.5, 0.1)){
  # Get data
  cohorts <- unique(cohort.names.short)
  point.ests <- as.numeric(data[csq, paste(cohorts, cancer, "AF", sep="_")])
  sample.sizes <- as.numeric(data[csq, paste(cohorts, cancer, "AN", sep="_")])
  meta.stats <- as.numeric(data[csq, paste("meta", cancer,
                                           c("AF", "AF_lower", "AF_upper"),
                                           sep="_")])

  # Get plot dimensions
  xlims <- c(0, 1 + (3 * x.space))
  if(is.null(ylims)){
    ylims <- c(0, max(point.ests, meta.stats))
  }
  if(proportional.bars){
    bar.xd <- c(0, cumsum(sample.sizes))/sum(sample.sizes)
  }else{
    bar.xd <- seq(0, 1, length=length(cohorts) + 1)
  }
  bar.xleft <- x.space + ((1 - diamond.ratio) * bar.xd[-length(bar.xd)])
  bar.xright <- x.space + ((1 - diamond.ratio) * bar.xd[-1])
  diamond.x <-  x.space + max(bar.xright) + (c(0, 0.5, 1, 0.5) * diamond.ratio)
  diamond.y <- meta.stats[c(1, 2, 1, 3)]

  # Prep plot area
  prep.plot.area(xlims=xlims, ylims=ylims, parmar=parmar)
  y.ax.at <- axTicks(2)
  if(length(y.ax.at) > 4){y.ax.at <- y.ax.at[c(TRUE, FALSE)]}
  if(label.y.axis){
    y.labels <- paste(100 * y.ax.at, "%", sep="")
    y.ax.col <- "black"
  }else{
    y.labels <- rep("", length(y.ax.at))
    y.ax.col <- "gray70"
  }
  if(title.y.axis){
    y.title <- cancer
  }else{
    y.title <- ""
  }
  if(y.axis){
    clean.axis(2, at=y.ax.at, tck=ax.tck, infinite=TRUE, labels=y.labels,
               title=y.title, title.line=1, cex.title=text.cex,
               cex.axis=text.cex, col.axis=y.ax.col)
  }
  mtext(title, side=3, font=2, cex=title.cex, line=title.line)

  # Add bars for individual cohorts
  rect(xleft=bar.xleft, xright=bar.xright, ybottom=0, ytop=point.ests,
       col=cancer.palettes[[cancer]][cohort.color.prefixes[cohorts]],
       border=cancer.palettes[[cancer]][cohort.color.prefixes[cohorts]])
  clean.axis(1, tck=0, labels=NA)

  # Add diamond for meta-analysis
  polygon(x=diamond.x, y=diamond.y, col=cancer.colors[cancer])
  segments(x0=diamond.x[1:2], x1=diamond.x[3:2],
           y0=meta.stats[1:2], y1=meta.stats[c(1, 3)],
           col="white")
  polygon(x=diamond.x, y=diamond.y, border="black", col=NA, xpd=T)
  text(x=diamond.x[1], y=diamond.y[4]-(0.035*diff(par("usr")[3:4])),
       pos=3, cex=text.cex, xpd=T,
       labels=paste(round(100 * meta.stats[1], 1), "%", sep=""))
}

# Simple legend for single cancer bar chart
simple.bar.legend <- function(cancer, text.cex=5/6,
                              parmar=c(0.1, 0.1, 0.1, 0.1)){
  cohorts <- unique(cohort.names.short)
  prep.plot.area(xlims=c(0, 4), ylim=c(4, 0), parmar=parmar)
  rect(xleft=rep(0.1, 3), xright=rep(0.9, 3), ybottom=(0:2)+0.1, ytop=(1:3)-0.1,
       col=cancer.palettes[[cancer]][cohort.color.prefixes[cohorts]])
  polygon(x=c(0.1, 0.5, 0.9, 0.5), y=3+c(0.5, 0.1, 0.5, 0.9),
          col=cancer.colors[cancer])
  segments(x0=c(0.1, 0.5), x1=c(0.9, 0.5), y0=3+c(0.5, 0.1), 3+c(0.5, 0.9), col="white")
  polygon(x=c(0.1, 0.5, 0.9, 0.5), y=3+c(0.5, 0.1, 0.5, 0.9))
  text(x=rep(0.5, 4), y=(0:3)+0.5, pos=4, labels=c(cohorts, "Meta."), cex=text.cex)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Summarize all KRAS somatic mutation data")
parser$add_argument("--stats", metavar=".tsv", type="character", required=TRUE,
                    help=paste("somatic mutation frequency .tsv generated by",
                               "gather_somatic_ras_data.py"))
parser$add_argument('--out-pdf', metavar='path', type="character", nargs=1,
                    help='output .pdf')
args <- parser$parse_args()

# # DEV
# args <- list("stats"="~/scratch/ras_somatic_variants.tsv.gz",
#              "out_pdf"="~/scratch/KRAS_summary_plot.pdf")

# Load somatic mutation frequency data
data <- load.somatic.freqs(args$stats)
data$csq.simple <- paste(data$ref, data$codon, data$alt, sep="")
data$max_freq <- apply(data[, grep("meta_.*_AF$", colnames(data))], 1, max)

# Restrict frequency data to top N non-CNA mutations
n.mut.highlights <- length(which(data$max_freq >= 0.01 & !(data$alt %in% c("DEL", "AMP"))))
freq.order <- order(data$max_freq, decreasing=TRUE)
highlight.muts <- head(setdiff(rownames(data)[freq.order],
                               c("AMP", "DEL", "any_mis")), n.mut.highlights)
top.data <- data[which(rownames(data) %in% highlight.muts), ]
top.data <- top.data[with(top.data, order(codon, -max_freq, alt)), ]
top.data <- rbind(top.data, data["any_mis", ])

# Set hard-coded plotting variables
n.codons <- 189
cancer.order <- c("PDAC", "CRAD", "LUAD")
ctype.ylims <- lapply(cancer.order, function(cancer){
  c(0, 1.25*max(top.data[which(rownames(top.data) != "any_mis"),
                         grep(paste(cancer, "AF", sep="_"), colnames(top.data))]))
})
names(ctype.ylims) <- cancer.order
protein.feat.colors <- c("P-Loop" = "#569DA8",
                         "GTP" = "#95D2DB",
                         "Switch I" = "#FAE641",
                         "Switch II" = "#FA9E41",
                         "Effector" = "#A85D5E",
                         "Hypervariable" = "#756A5E")

# Plot main figure
pdf(args$out_pdf, width=180/25.4, height=4.2)
layout(matrix(c(rep(1, n.mut.highlights+2),
                rep(2, n.mut.highlights+2),
                rep(3, n.mut.highlights), 4, 4,
                rep(5, n.mut.highlights), 6, 6,
                rep(7, n.mut.highlights), 10, 10,
                rep(8, n.mut.highlights), 10, 10,
                rep(9, n.mut.highlights), 10, 10,
                (1:(3*(n.mut.highlights+2)))+10),
              nrow=10, byrow=T),
       heights=c(0.4, 1, 0.45, 0.5, 0.15, 0.15, 0.15, 1.2, 0.75, 0.75),
       widths=c(1.4, rep(1, n.mut.highlights-1), 1.3, 1.1))
kras.idio(parmar=c(0.5, 2.5, 0.5, 0.5))
kras.gene(parmar=c(0.2, 2.5, 1.25, 0.5))
kras.protein(label.domains=FALSE, parmar=c(0.5, 3.2, 0.5, 0.5))
protein.legend()
codon.heat(counts=get.codon.mut.counts(data),
           pal=colorRampPalette(stage.colors)(6),
           left.label="Residues",
           parmar=c(0.5, 3.2, 1, 0.5))
axis(3, at=c(12, 60.5, 145.5), tick=F, line=-1, cex.axis=5/6,
     labels=c("G12-G13", "Q61", "A146"))
codon.mut.heat.legend(parmar=c(0.5, 0.5, 0.5, 0.5))
codon.heat(counts=get.codon.mut.freq(data, "PDAC"),
           pal=viridis(51), left.label="PDAC",
           number.residues=FALSE, parmar=c(0.1, 3.2, 0.1, 0.5))
codon.heat(counts=get.codon.mut.freq(data, "CRAD"),
           pal=viridis(51), left.label="CRAD",
           number.residues=FALSE, parmar=c(0.1, 3.2, 0.1, 0.5))
codon.heat(counts=get.codon.mut.freq(data, "LUAD"),
           pal=viridis(51), left.label="LUAD",
           number.residues=FALSE, parmar=c(0.1, 3.2, 0.1, 0.5))
codon.mut.freq.heat.legend(parmar=c(0.8, 1.5, 0.2, 2.5))
for(k in 1:length(cancer.order)){
  for(i in 1:n.mut.highlights){
    parmar <- if(i==1){c(0.5, 3, 0.35, 0.5)}else{c(0.5, 1, 0.35, 0.5)}
    if(k == 1){parmar[3] <- 3.25}
    plot.freq.meta.single(top.data, rownames(top.data)[i], cancer=cancer.order[k],
                          ylims=ctype.ylims[[k]], label.y.axis=(i == 1), title.y.axis=(i == 1),
                          title=if(k==1){top.data[i, "csq.simple"]}else{NULL},
                          parmar=parmar)
  }
  plot.freq.meta.single(top.data, "any_mis", cancer=cancer.order[k],
                        ylims=c(0, 1), title=if(k==1){"Any Mis."}else{NULL},
                        title.y.axis=FALSE,
                        parmar=if(k==1){c(0.5, 2.5, 3.25, 0.5)}else{c(0.5, 2.5, 0.35, 0.5)})
  simple.bar.legend(cancer.order[k],
                    parmar=if(k==1){c(0.75, 0.5, 3.25, 0.5)}else{c(0.75, 0.5, 0.35, 0.5)})
}
dev.off()

