#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Fit multivariate regression model to estimate marginal effects of FGFR4 common
# and rare germline variants on somatic KRAS tier 1 mutation status

# Also compiles a table of counts for KRAS alleles by FGFR4 genotype


#########
# Setup #
#########
# Load necessary libraries and constants
require(argparse, quietly=TRUE)
require(OncoModR, quietly=TRUE)
OncoModR::load.constants("all")

# Declare FGFR4 protein domain coordinates in amino acid space
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
kras.aa.abbrevs <- c("Gly12Ala" = "G12A",
                     "Gly12Asp" = "G12D",
                     "Gly12Cys" = "G12C",
                     "Gly12Ser" = "G12R",
                     "Gly12Val" = "G12V",
                     "Gly13Asp" = "G13D")


##################
# Data Functions #
##################
# Build a map of KRAS aa to nt changes
make.aa.map <- function(set_map_tsv, somatic.vids){
  set.map <- read.table(set_map_tsv, header=F, sep="\t")
  colnames(set.map) <- c("csq", "vids")
  set.map$vids <- sapply(set.map$vids, function(mstr){unlist(strsplit(mstr, split=","))})
  ovr.idx <- which(sapply(set.map$vids, function(vids){length(intersect(vids, somatic.vids)) > 0}))
  aa.idx <- grep("_p.", set.map$csq, fixed=T)
  set.map <- set.map[intersect(ovr.idx, aa.idx), ]
  set.map$csq <- sapply(set.map$csq, function(s){unlist(strsplit(s, split="_p.", fixed=T))[2]})
  return(set.map[order(kras.aa.abbrevs[set.map$csq]), ])
}

# Compile a table of KRAS allele counts vs. cohort
count.kras.alleles.by.cohort <- function(somatic.ad, cohort.names, kras.aa.map,
                                         keep.samples=NULL){
  res <- as.data.frame(do.call("cbind", lapply(1:nrow(kras.aa.map), function(aidx){
    sapply(somatic.ad, function(df){
      if(is.null(keep.samples)){
        keep.samples <- colnames(df)
      }
      ad.vals <- apply(df[intersect(rownames(df), unlist(kras.aa.map[aidx, 2])),
         intersect(colnames(df), keep.samples)], 2, sum, na.rm=T)
      length(which(ad.vals > 0))
    })
  })))
  colnames(res) <- kras.aa.abbrevs[kras.aa.map[, 1]]
  rownames(res) <- cohort.names
  return(res)
}


######################
# Plotting Functions #
######################
# Plot coefficients from multivariate regression model
plot.multivariate.coeffs <- function(fit.df){
  coeffs.to.plot <- c("H3", "V10I", "G388R", "nonsyn_key_domains",
                      "nonsyn_cytoplasmic", "nonsyn_other", "rare_synonymous")
  plot.df <- fit.df[intersect(coeffs.to.plot, rownames(fit.df)), ]
  plot.df$lower <- plot.df$Estimate + (qnorm(0.025) * plot.df$`Std. Error`)
  plot.df$upper <- plot.df$Estimate + (qnorm(0.975) * plot.df$`Std. Error`)
  plot.df <- plot.df[order(-plot.df[, 4]), ]
  plot.df[, c("Estimate", "lower", "upper")] <- log2(exp(plot.df[, c("Estimate", "lower", "upper")]))
  plot.labels <- c("G388R" = "G388R",
                   "nonsyn_key_domains" = "Rare nonsynonymous\nin transmem. domain\nor at disulfide bonds",
                   # "H3" = "Haplotype 3\n(includes P136L)",
                   "V10I" = "V10I",
                   "nonsyn_cytoplasmic" = "Rare nonsynonymous\nin cytoplasmic domain",
                   "nonsyn_other" = "Other rare\nnonsynonymous",
                   "rare_synonymous" = "Rare\nsynonymous")
  plot.colors <- c("G388R" = as.character(csq.colors["missense"]),
                   "nonsyn_key_domains" = as.character(csq.colors["missense"]),
                   # "H3" = as.character(csq.colors["missense"]),
                   "V10I" = as.character(csq.colors["missense"]),
                   "nonsyn_cytoplasmic" = as.character(csq.colors["missense"]),
                   "nonsyn_other" = as.character(csq.colors["missense"]),
                   "rare_synonymous" = as.character(csq.colors["synonymous"]))
  plot.colors.shaded <- sapply(plot.df$Variable, function(v){
    alpha <- if(plot.df[v, "Pr(>|z|)"] <= 0.05){1}else{0.3}
    adjustcolor(as.character(plot.colors[v]), alpha=alpha)
  })
  plot.label.colors <- c("TRUE"="black", "FALSE"="gray75")[as.character(plot.df$`Pr(>|z|)` <= 0.05)]
  prep.plot.area(xlims=c(min(c(log2(0.75), 1.5*min(plot.df$Estimate))), 1.5*max(plot.df$Estimate)),
                 ylims=c(nrow(plot.df), 0), parmar=c(0.25, 7, 3, 4))
  abline(v=0)
  rect(xleft=0, xright=plot.df$Estimate,
       ybottom=(1:nrow(plot.df))-0.9, ytop=(1:nrow(plot.df))-0.1,
       col=plot.colors.shaded, border=plot.colors[plot.df$Variable])
  segments(x0=plot.df$lower, x1=plot.df$upper,
           y0=(1:nrow(plot.df))-0.5, y1=(1:nrow(plot.df))-0.5,
           col=plot.label.colors)
  sapply(1:nrow(plot.df), function(y){
    axis(2, at=y-0.5, labels=plot.labels[plot.df$Variable[y]],
         tick=F, line=-1, cex.axis=5/6, col.axis=plot.label.colors[y],
         las=2)
  })
  sapply(1:nrow(plot.df), function(y){
    axis(4, at=y-0.5, labels=RLCtools::format.pval(plot.df$`Pr(>|z|)`[y]),
         tick=F, line=-1, cex.axis=5/6, col.axis=plot.label.colors[y],
         las=2)
  })
  l2.at <- log2(c(0.75, 1, 1.25, 1.5, 2, 3, 4, 5))
  l2.labels <- c(0.75, 1, 1.25, 1.5, 2, 3, 4, 5)
  clean.axis(3, at=l2.at, infinite=TRUE, labels=paste(l2.labels, "x", sep=""),
             title=bquote("Odds of somatic" ~ italic("KRAS") ~ "mutation per germline" ~ italic("FGFR4") ~ "allele"))
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Custom multivariate regression model",
                                           "for germline FGFR4 alleles vs. somatic",
                                           "KRAS tier 1 mutations"))
parser$add_argument("--sample-metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv. One per cohort, in order.",
                    required=TRUE, action="append")
parser$add_argument("--somatic-ad", metavar=".tsv", type="character",
                    help="Somatic allele dosage matrix. One per cohort, in order.",
                    required=TRUE, action="append")
parser$add_argument("--germline-ad", metavar=".tsv", type="character",
                    help="Germline allele dosage matrix. One per cohort, in order.",
                    required=TRUE, action="append")
parser$add_argument("--name", metavar="string", type="character", action="append",
                    help=paste("names for each set of cohort inputs, provided",
                               "in the same order"))
parser$add_argument("--somatic-variant-ids", metavar=".txt", type="character",
                    help="List of somatic mutation IDs to consider tier 1", required=TRUE)
parser$add_argument("--somatic-variant-sets", metavar=".tsv", type="character",
                    help="Map of --somatic-variant-ids to somatic endpoints",
                    required=TRUE)
parser$add_argument("--germline-haplotype-snps", metavar=".tsv", type="character",
                    help="Two-column .tsv mapping variant IDs to haplotypes",
                    required=TRUE)
parser$add_argument("--coding-variant-map", metavar=".tsv", type="character",
                    help="Two-column .tsv of FGFR4 protein variants to germline variant id(s)", required=TRUE)
parser$add_argument("--eligible-controls", metavar="path", type="character",
                    help=paste("path to one or more lists of samples eligible to",
                               "be treated as controls for association testing.",
                               "If no file is provided, all samples are eligible."),
                    action="append")
parser$add_argument("--outdir", metavar="path", type="character",
                    help="output directory for results", default="./")
parser$add_argument("--cancer-type", metavar="character",
                    help=paste("Subset to samples from this cancer type",
                               "[default: use all samples]"))
args <- parser$parse_args()

# # DEV:
# args <- list("sample_metadata" = c("~/scratch/TCGA.ALL.sample_metadata.tsv.gz",
#                                    "~/scratch/PROFILE.ALL.sample_metadata.tsv.gz",
#                                    "~/scratch/HMF.ALL.sample_metadata.tsv.gz"),
#              "somatic_ad" = c("~/scratch/TCGA.somatic_variants.dosage.tsv.gz",
#                               "~/scratch/PROFILE.somatic_variants.dosage.tsv.gz",
#                               "~/scratch/HMF.somatic_variants.dosage.tsv.gz"),
#              "germline_ad" = c("~/scratch/TCGA.FGFR4.dosage.tsv.gz",
#                                "~/scratch/PROFILE.FGFR4.dosage.tsv.gz",
#                                "~/scratch/HMF.FGFR4.dosage.tsv.gz"),
#              "name" = c("TCGA", "PROFILE", "HMF"),
#              "somatic_variant_ids" = "~/scratch/KRAS_tier_1.somatic.variant_ids.list",
#              "somatic_variant_sets" = "~/scratch/CRAD.KRAS.somatic_endpoints.tsv",
#              "germline_haplotype_snps" = "~/scratch/FGFR4.common_variants.haplotype_assignment.tsv",
#              "coding_variant_map" = "~/scratch/FGFR4.coding_variants.sets.tsv",
#              "eligible_controls" = c("~/scratch/TCGA.ALL.eligible_controls.list",
#                                      "~/scratch/PROFILE.ALL.eligible_controls.list",
#                                      "~/scratch/HMF.ALL.eligible_controls.list"),
#              "outdir" = "~/scratch/",
#              "cancer_type" = "CRAD")

# Sanity check to make sure all cohorts have the same number of inputs
if(is.null(args$name)){
  args$name <- paste("cohort", 1:length(args$stats), sep="")
}
cohort.names <- args$name
n.cohorts <- length(cohort.names)
if(length(args$sample_metadata) != n.cohorts |
   length(args$somatic_ad) != n.cohorts |
   length(args$germline_ad) != n.cohorts){
  stop("Different number of inputs for --sample-metadata, --somatic-ad, and --germline-ad. Exiting.")
}

# Load patient metadata, subset to cancer type of interest (if optioned), and merge
meta.list <- lapply(args$sample_metadata, load.patient.metadata, deduplicate=TRUE)
meta <- merge.patient.metadata(meta.list, args$cancer_type, args$name)
samples.w.pheno <- rownames(meta)

# Impute cohort-specific PCs in all other cohorts as mean
# Hopefully, this will allow cohort-specific ancestry weights
# to be fit by the association models
meta <- impute.missing.values(meta, fill.missing="mean")

# Load lists of eligible control samples, if optioned
all.elig.controls <- unique(unlist(sapply(args$eligible_controls, function(path){
  as.character(read.table(path)[, 1])
})))
elig.controls <- intersect(samples.w.pheno, all.elig.controls)

# Load germline and somatic variant lists
somatic.vids <- sort(unique(read.table(args$somatic_variant_ids, header=F)[, 1]))
kras.aa.map <- make.aa.map(args$somatic_variant_sets, somatic.vids)
hap.snps <- read.table(args$germline_haplotype_snps, header=T, sep="\t", comment.char="")
colnames(hap.snps) <- c("vid", "haplotype")
coding.variants <- load.variant.sets(args$coding_variant_map)
coding.variants$variant_ids <- sapply(coding.variants$variant_ids, function(vstr){unlist(strsplit(vstr, split=","))})
coding.variants$csq <- "missense"
coding.variants$csq[grep("=$", coding.variants$set_id)] <- "synonymous"
coding.variants$csq[grep("Ter$|Ter[0-9]+$", coding.variants$set_id)] <- "lof"
coding.variants$aa <- sapply(coding.variants$set_id, function(gstr){
  as.numeric(sub("=", "", unlist(strsplit(unlist(strsplit(gstr, split="_p.", fixed=T))[2], split="[A-z]+"))[2]))
})

# Load germline and somatic allele depth matrixes
germline.ad <- lapply(args$germline_ad, load.ad.matrix, sample.subset=samples.w.pheno)
somatic.ad <- lapply(args$somatic_ad, load.ad.matrix, sample.subset=samples.w.pheno,
                     variant.subset=somatic.vids)

# # Make copies of germline AD matrixes with "simple" variant IDs
# # This is required for haplotype dosage estimation (below)
# germline.ad.simpleIDs <- lapply(germline.ad, simplify.vids)
#
# # Estimate haplotype dosage for each sample based on available data
# haplotypes <- sort(unique(hap.snps$haplotype))
# hap.dosage <- do.call("cbind", lapply(haplotypes, function(hid){
#   marker.vids <- hap.snps$vid[which(hap.snps$haplotype == hid)]
#   query.ad.matrix(germline.ad.simpleIDs, vids=marker.vids, action="mean")
# }))
# colnames(hap.dosage) <- haplotypes
# hap.dosage <- impute.missing.values(hap.dosage, fill.columns=haplotypes)
# hap.dosage <- as.data.frame(t(apply(hap.dosage, 1, function(h.est){
#   if(sum(h.est) < 2){
#     h.est <- h.est + ((2 - sum(h.est)) / length(haplotypes))
#   }
#   2 * h.est / sum(h.est, na.rm=T)
# })))
# rm(germline.ad.simpleIDs)

# Estimate allele frequency for all coding variants
coding.afs <- lapply(coding.variants$variant_ids, function(vlist){
  v <- mean(query.ad.matrix(germline.ad, vids=unlist(vlist), action="sum"), na.rm=T) / 2
  if(!is.na(v)){v}else{0}
})
coding.variants$AF <- coding.afs

# Get germline allele counts for various sets of FGFR4 alleles of interest
h3.idx <- which(coding.variants$set_id == "ENST00000292408_p.Pro136Leu")
h3.ac <- query.ad.matrix(germline.ad,
                         vids=unlist(coding.variants$variant_ids[h3.idx]),
                         action="sum", missing.vid.fill=0)
g388r.idx <- which(coding.variants$set_id == "ENST00000292408_p.Gly388Arg")
g388r.ac <- query.ad.matrix(germline.ad,
                           vids=unlist(coding.variants$variant_ids[g388r.idx]),
                           action="sum", missing.vid.fill=0)
v10i.idx <- which(coding.variants$set_id == "ENST00000292408_p.Val10Ile")
v10i.ac <- query.ad.matrix(germline.ad,
                           vids=unlist(coding.variants$variant_ids[v10i.idx]),
                           action="sum", missing.vid.fill=0)
# rare_nonsyn.idxs <- which(coding.variants$csq != "synonymous"
#                               & coding.variants$AF < 0.01)
# rare_nonsyn.vids <- unlist(coding.variants$variant_ids[rare_nonsyn.idxs])
# rare_nonsyn.ac <- query.ad.matrix(germline.ad, vids=rare_nonsyn.vids,
#                                       action="sum", missing.vid.fill=0)
key_domains.idxs <- which(coding.variants$csq != "synonymous"
                          & coding.variants$AF < 0.01
                          & (sapply(coding.variants$aa, is.inside, interval=protein.regions[[2]])
                             | sapply(coding.variants$aa, is.inside, interval=disulfide.bonds[[2]])
                             | sapply(coding.variants$aa, is.inside, interval=disulfide.bonds[[3]])))
key_domains.vids <- unlist(coding.variants[key_domains.idxs, "variant_ids"])
key_domains.ac <- query.ad.matrix(germline.ad, vids=key_domains.vids,
                                  action="sum", missing.vid.fill=0)
cytoplasmic.idxs <- which(coding.variants$csq != "synonymous"
                          & coding.variants$AF < 0.01
                          & sapply(coding.variants$aa, is.inside, interval=protein.regions[[3]]))
cytoplasmic.vids <- unlist(coding.variants$variant_ids[cytoplasmic.idxs])
cytoplasmic.ac <- query.ad.matrix(germline.ad, vids=cytoplasmic.vids,
                                  action="sum", missing.vid.fill=0)
nonsyn_other.idxs <- setdiff(which(coding.variants$csq != "synonymous"
                                   & coding.variants$AF < 0.01),
                             c(key_domains.idxs, cytoplasmic.idxs))
nonsyn_other.vids <- unlist(coding.variants$variant_ids[nonsyn_other.idxs])
nonsyn_other.ac <- query.ad.matrix(germline.ad, vids=nonsyn_other.vids,
                                   action="sum", missing.vid.fill=0)
rare_synonymous.idxs <- which(coding.variants$csq == "synonymous"
                                   & coding.variants$AF < 0.01)
rare_synonymous.vids <- unlist(coding.variants$variant_ids[rare_synonymous.idxs])
rare_synonymous.ac <- query.ad.matrix(germline.ad, vids=rare_synonymous.vids,
                                   action="sum", missing.vid.fill=0)

# Get somatic mutation status
Y.vals <- query.ad.matrix(somatic.ad, somatic.vids, elig.controls=elig.controls, action="any", missing.vid.fill=0)
Y.vals <- Y.vals[which(!is.na(Y.vals))]

# Build test dataframe
meta.samples <- rownames(meta)
Y.samples <- names(Y.vals)
# haplotype.samples <- rownames(hap.dosage)
# final.samples <- intersect(intersect(meta.samples, Y.samples), haplotype.samples)
final.samples <- intersect(meta.samples, Y.samples)
test.df <- meta[final.samples, c("SEX", "AGE_AT_DIAGNOSIS", "TUMOR_PURITY",
                      colnames(meta)[grep("^cohort\\.", colnames(meta))],
                      colnames(meta)[grep("^PC[1-3]\\.", colnames(meta))])]
test.df$SEX <- as.numeric(test.df$SEX == "MALE")
test.df$Y <- Y.vals[final.samples]
# test.df$H3 <- hap.dosage[final.samples, "H3"]
# test.df$H3 <- h3.ac[final.samples]
# test.df$nonH3 <- 2 - h3.ac[final.samples]
test.df$V10I <- v10i.ac[final.samples]
test.df$G388R <- g388r.ac[final.samples]
# test.df$rare_nonsyn <- rare_nonsyn.ac[final.samples]
test.df$nonsyn_key_domains <- key_domains.ac[final.samples]
test.df$nonsyn_cytoplasmic <- cytoplasmic.ac[final.samples]
test.df$nonsyn_other <- nonsyn_other.ac[final.samples]
test.df$rare_synonymous <- rare_synonymous.ac[final.samples]
# test.df$heterozygous_haplotypes <- test.df$V10I == 1 | test.df$G388R == 1 | round(test.df$H3) == 1

# Standard normalize all non-integer/non-factor covariates
no.scale.covariates <- c("Y", "SEX",
                         colnames(test.df)[grep("^cohort\\.", colnames(test.df))],
                         "H3", "nonH3", "V10I", "G388R",
                         colnames(test.df)[grep("syn", colnames(test.df))])
scale.cidxs <- which(!(colnames(test.df) %in% no.scale.covariates))
test.df[, scale.cidxs] <- apply(test.df[, scale.cidxs], 2, scale)

# Fit regression
fit <- glm(Y ~ ., data=test.df, family="binomial")
fit.df <- as.data.frame(summary(fit)$coefficients)
fit.df$Variable <- rownames(fit.df)
fit.df <- fit.df[, c("Variable", setdiff(colnames(fit.df), "Variable"))]

# Write summary of regression coefficients
write.table(fit.df,
            paste(args$outdir, "FGFR4.multivariate_regression_results.tsv", sep="/"),
            col.names=T, row.names=F, quote=F, sep="\t")

# Plot selected coefficients
pdf(paste(args$outdir, "FGFR4.multivariate_association.coefficients.pdf", sep="/"),
    height=4.2, width=6)
plot.multivariate.coeffs(fit.df)
dev.off()

# Cross-tabs of EUR-only for two significant germline categories
eur.ids <- rownames(meta)[which(meta$POPULATION == "EUR")]
lapply(1:length(meta.list), function(i){
  cohort.ids <- meta.list[[i]][, 1]
  sub.ids <- intersect(rownames(test.df), intersect(eur.ids, cohort.ids))
  print(table(test.df[sub.ids, c("Y", "G388R")]))
})

# Cross-tabs of KRAS allele rates by FGFR genotype & cohort
cat("\nRare missense carriers in key FGFR4 domains:\n")
key_domains.sids <- rownames(test.df)[which(test.df$nonsyn_key_domains > 0)]
key_domains.aa_rate <- count.kras.alleles.by.cohort(somatic.ad, cohort.names,
                                                    kras.aa.map,
                                                    keep.samples=key_domains.sids)
print(key_domains.aa_rate)
cat("\nFGFR4 G388R carriers:\n")
g388r.sids <- rownames(test.df)[which(test.df$nonsyn_key_domains < 1
                                      & test.df$G388R > 0)]
g388r.aa_rate <- count.kras.alleles.by.cohort(somatic.ad, cohort.names,
                                              kras.aa.map, keep.samples=g388r.sids)
print(g388r.aa_rate)
cat("\nAll other colorectal samples:\n")
other.aa_rate <- count.kras.alleles.by.cohort(somatic.ad, cohort.names, kras.aa.map,
                                   keep.samples=setdiff(rownames(test.df), c(key_domains.sids, g388r.sids)))
print(other.aa_rate)


# # Statistical test of key_domains vs. everything else
# chisq.test(cbind(apply(key_domains.aa_rate, 2, sum),
# apply(rbind(g388r.aa_rate, other.aa_rate), 2, sum)))

