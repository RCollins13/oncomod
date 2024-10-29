#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Fit multivariate regression model to attempt replication of FGFR4-KRAS
# association based on data from Gurjao et al., Cancer Disc., 2021


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
tier1.muts <- c("G12D", "G12V", "G13D", "G12R", "G12A", "G12C")
control.exclusion.list <- c("G12A", "G12C", "G12D", "G12R", "G12S", "G12T", "G12V",
                            "Q61H", "Q61K", "Q61L", "Q61R", "A146P", "A146T",
                            "A146V", "G10dup", "G13C", "G13D", "G13dup", "G13R",
                            "G13V", "L19F", "Q22K", "D33E", "A59T", "K117N", "K117R")


##################
# Data Functions #
##################
# Load sample metadata and restrict to samples for inclusion in model
load.gurjao.meta <- function(meta.in, pcs.in){
  # Load metadata
  meta <- read.table(meta.in, header=T, sep="\t", comment.char="")

  # Reformat/refactor variables for inclusion as covariates
  meta$Cohort <- relevel(as.factor(meta$Cohort), ref="NHS")
  meta$Age <- as.numeric(meta$Age)
  meta <- impute.missing.values(meta, fill.columns=c("Age", "Sex"))
  meta$advanced_disease <- meta$Stage %in% c("III", "IV")

  # Load & append PCs
  pcs <- read.table(pcs.in, sep="\t", comment.char="", header=T)
  rownames(pcs) <- sapply(strsplit(pcs[, 1], split="-"), function(l){l[1]})
  pcs[, 1] <- NULL
  meta <- merge(meta, pcs, by.x="ID", by.y="row.names", all.x=T, all.y=F, sort=F)

  # Return only required information for MSS tumors
  meta[which(meta$Microsatellite.status == "MSS"),
       c("ID", "Cohort", "Age", "Sex", "advanced_disease",
         colnames(meta)[grep("^PC[1-9]", colnames(meta))])]
}

# Add KRAS mutation status to sample metadata and only retain eligible controls or Tier 1 mutants
add.kras.status <- function(meta, tsv.in){
  # Load & clean KRAS mutations
  kras.df <- read.table(tsv.in, header=T, sep="\t", comment.char="")
  kras.df$ID <- do.call("rbind", strsplit(kras.df$pair_id, split="-", fixed=T))[, 1]
  kras.df$aa_csq <- sub("^p.", "", kras.df$aa_csq)

  # Get list of samples with tier 1 mutation or other mutation
  # Assume all other samples are WT (elig controls)
  tier1.ids <- unique(kras.df$ID[which(kras.df$aa_csq %in% tier1.muts)])
  noncontrol.ids <- unique(kras.df$ID[which(kras.df$aa_csq %in% control.exclusion.list)])
  wt.ids <- setdiff(meta$ID, noncontrol.ids)

  # Annotate sample metadata with KRAS mutation status
  # Only keep WT or Tier 1 samples
  meta$tier1 <- meta$ID %in% tier1.ids
  meta$wt <- meta$ID %in% wt.ids
  meta[which(meta$tier1 | meta$wt), ]
}

# Load & parse table of FGFR4 variant consequences
load.csq.table <- function(tsv.in){
  csq.df <- read.table(tsv.in, header=T, sep="\t", comment.char="")
  colnames(csq.df)[1] <- "VID"
  csq.df[, c("aa_ref", "aa_alt")] <- do.call("rbind", strsplit(csq.df$aa_change, split="/"))
  csq.df$csq <- sub("_variant$", "", csq.df$csq)
  csq.df$csq[which(csq.df$csq %in% c("frameshift", "splice_donor"))] <- "lof"
  csq.df[which(csq.df$csq %in% names(csq.names.short)), ]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Multivariate regression of germline ",
                                           "FGFR4 alleles vs. somatic",
                                           "KRAS tier 1 mutations in data from ",
                                           "Gurjao et al., 2021"))
parser$add_argument("--sample-metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--kras-status", metavar=".tsv", type="character",
                    help="Table of reported somatic KRAS mutations",
                    required=TRUE)
parser$add_argument("--germline-ad", metavar=".tsv", type="character",
                    help="Germline allele dosage matrix", required=TRUE)
parser$add_argument("--germline-variant-csq", metavar=".tsv", type="character",
                    help="Table of consequences for all germline variants",
                    required=TRUE)
parser$add_argument("--pcs", metavar=".tsv", type="character", required=TRUE,
                    help=paste("Principal components for all samples computed",
                               "from common SNPs"))
parser$add_argument("--outdir", metavar="path", type="character",
                    help="output directory for results", default="./")
args <- parser$parse_args()

# # DEV:
# args <- list("sample_metadata" = "~/Desktop/Collins/VanAllen/RAS_projects/RAS_modifiers/FGFR4_replication/giannakis/gurjao_2021.table_s1.demographics.tsv",
#              "kras_status" = "~/Desktop/Collins/VanAllen/RAS_projects/RAS_modifiers/FGFR4_replication/giannakis/gurjao_2021.kras_mutations.tsv.gz",
#              "germline_ad" = "~/Desktop/Collins/VanAllen/RAS_projects/RAS_modifiers/FGFR4_replication/giannakis/gurjao_2021.FGFR4.vep.clean.gq10.dosage.tsv.gz",
#              "germline_variant_csq" = "~/Desktop/Collins/VanAllen/RAS_projects/RAS_modifiers/FGFR4_replication/giannakis/gurjao_2021.FGFR4.gq10.variant_csqs.tsv",
#              "pcs" = "~/Desktop/Collins/VanAllen/RAS_projects/RAS_modifiers/FGFR4_replication/giannakis/gurjao_2021.common_SNPs.eigenvec",
#              "outdir" = "~/scratch/")

# Load sample metadata and filter to samples for consideration in analysis
meta <- load.gurjao.meta(args$sample_metadata, args$pcs)

# Dev code to quickly estimate European samples (those with all PC Z-scores within ~ [-1, 1])
# eur.ids <- rownames(meta)[which(apply(apply(meta[, grep("^PC", colnames(meta))], 2, scale), 1, function(v){all(v < 1 & v > -1)}))]
# meta <- meta[eur.ids, ]

# Add KRAS mutation status to sample metadata and only retain eligible controls or Tier 1 mutants
meta <- add.kras.status(meta, args$kras_status)

# Load germline variant consequences
csq.df <- load.csq.table(args$germline_variant_csq)

# Load germline allele depth matrix and correct sample IDs
germline.ad <- load.ad.matrix(args$germline_ad)
colnames(germline.ad) <- sapply(strsplit(colnames(germline.ad), split="-", fixed=T), function(l){l[1]})
keep.samples <- intersect(meta$ID, colnames(germline.ad))
cat(paste("Retained", prettyNum(length(keep.samples), big.mark=","), "patients for replication\n"))
germline.ad <- germline.ad[, keep.samples]
rownames(meta) <- meta$ID
meta$ID <- NULL
meta <- meta[keep.samples, ]

# Get germline allele counts for various sets of FGFR4 alleles of interest
all_nonsyn.idxs <- which(csq.df$csq != "synonymous")
all_nonsyn.vids <- unlist(csq.df$VID[all_nonsyn.idxs])
all_nonsyn.ac <- query.ad.matrix(germline.ad, vids=all_nonsyn.vids, action="sum")
g388r.idx <- which(csq.df$aa_number == 388 & csq.df$aa_change == "G/R")
g388r.ac <- query.ad.matrix(germline.ad, vids=unlist(csq.df$VID[g388r.idx]), action="sum")
v10i.idx <- which(csq.df$aa_number == 10 & csq.df$aa_change == "V/I")
v10i.ac <- query.ad.matrix(germline.ad, vids=unlist(csq.df$VID[v10i.idx]), action="sum")
rare_nonsyn.idxs <- which(csq.df$csq != "synonymous" & csq.df$AF < 0.01)
rare_nonsyn.vids <- unlist(csq.df$VID[rare_nonsyn.idxs])
rare_nonsyn.ac <- query.ad.matrix(germline.ad, vids=rare_nonsyn.vids,
                                  action="sum", missing.vid.fill=0)
rare_synonymous.idxs <- which(csq.df$csq == "synonymous"
                              & csq.df$AF < 0.01)
rare_synonymous.vids <- unlist(csq.df$VID[rare_synonymous.idxs])
rare_synonymous.ac <- query.ad.matrix(germline.ad, vids=rare_synonymous.vids,
                                      action="sum", missing.vid.fill=0, )

# Get somatic mutation status
Y.vals <- as.numeric(meta$tier1)
names(Y.vals) <- rownames(meta)

# Build test dataframe
meta.samples <- rownames(meta)
Y.samples <- names(Y.vals)
final.samples <- intersect(meta.samples, Y.samples)
test.df <- meta[final.samples,
                c("Cohort", "Age", "Sex", "advanced_disease",
                  colnames(meta)[grep("^PC[1-4]$", colnames(meta))])]
test.df$Y <- Y.vals[final.samples]
test.df$V10I <- v10i.ac[final.samples]
test.df$G388R <- g388r.ac[final.samples]
test.df$rare_nonsyn <- rare_nonsyn.ac[final.samples]
test.df$rare_synonymous <- rare_synonymous.ac[final.samples]
test.df$Age <- scale(test.df$Age)
pc.idxs <- grep("^PC[1-9]", colnames(test.df))
test.df[, pc.idxs] <- apply(test.df[, pc.idxs], 2, scale)
test.df <- impute.missing.values(test.df, fill.columns=c("rare_nonsyn", "rare_synonymous"))
# test.df <- test.df[complete.cases(test.df), ]

# Fit regression
fit <- glm(Y ~ ., data=test.df, family="binomial")
fit.df <- as.data.frame(summary(fit)$coefficients)
fit.df$Variable <- rownames(fit.df)
fit.df <- fit.df[, c("Variable", setdiff(colnames(fit.df), "Variable"))]

# Write summary of regression coefficients
write.table(fit.df,
            paste(args$outdir, "gurjao_2021.FGFR4.multivariate_regression_results.tsv", sep="/"),
            col.names=T, row.names=F, quote=F, sep="\t")
