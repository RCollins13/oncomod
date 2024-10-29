#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2024-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Perform survival analysis of cancer patients conditional on
# FGFR4 germline and KRAS somatic configurations


#########
# Setup #
#########
# Load necessary libraries and constants
require(argparse, quietly=TRUE)
require(OncoModR, quietly=TRUE)
require(survival, quietly=TRUE)
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

# Fit a Cox regression to the provided set of data
fit.cox <- function(surv.df, other.terms=NULL){
  keep.cols <- c("DAYS_SURVIVED", "DEAD", "SEX", "AGE_AT_DIAGNOSIS")
  surv.df <- surv.df[, intersect(unique(c(keep.cols, other.terms)), colnames(surv.df))]
  coxph(Surv(DAYS_SURVIVED, DEAD) ~ ., data=surv.df)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Perform survival analysis of cancer",
                                           "patients conditional on FGFR4 germline",
                                           "and KRAS somatic configurations"))
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
#
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

# Load & merge patient metadata
# Subset to European patients with survival data for cancer type of interest (if optioned)
# Impute missing values on a cohort-specific basis after filtering
meta.list <- lapply(args$sample_metadata, load.patient.metadata, deduplicate=TRUE)
meta.list <- lapply(meta.list, function(df){
  impute.missing.values(df[which(!is.na(df$DAYS_SURVIVED) & df$POPULATION == "EUR"), ])
  })
meta <- merge.patient.metadata(meta.list, args$cancer_type, args$name,
                               impute.missing=FALSE)
samples.w.surv <- rownames(meta)

# Load lists of eligible control samples, if optioned
all.elig.controls <- unique(unlist(sapply(args$eligible_controls, function(path){
  as.character(read.table(path)[, 1])
})))
elig.controls <- intersect(samples.w.surv, all.elig.controls)

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
germline.ad <- lapply(args$germline_ad, load.ad.matrix, sample.subset=samples.w.surv)
somatic.ad <- lapply(args$somatic_ad, load.ad.matrix, sample.subset=samples.w.surv,
                     variant.subset=somatic.vids)

# Estimate allele frequency for all coding variants
coding.afs <- sapply(coding.variants$variant_ids, function(vlist){
  mean(query.ad.matrix(germline.ad, vids=unlist(vlist), action="mean"), na.rm=T)
})
coding.variants$AF <- coding.afs

# Get germline allele counts for various sets of FGFR4 alleles of interest
g388r.idx <- which(coding.variants$set_id == "ENST00000292408_p.Gly388Arg")
g388r.ac <- query.ad.matrix(germline.ad,
                            vids=unlist(coding.variants$variant_ids[g388r.idx]),
                            action="sum", missing.vid.fill=0)
v10i.idx <- which(coding.variants$set_id == "ENST00000292408_p.Val10Ile")
v10i.ac <- query.ad.matrix(germline.ad,
                           vids=unlist(coding.variants$variant_ids[v10i.idx]),
                           action="sum", missing.vid.fill=0)
rare_nonsyn.idxs <- which(coding.variants$csq != "synonymous"
                              & coding.variants$AF < 0.01)
rare_nonsyn.vids <- unlist(coding.variants$variant_ids[rare_nonsyn.idxs])
rare_nonsyn.ac <- query.ad.matrix(germline.ad, vids=rare_nonsyn.vids,
                                      action="sum", missing.vid.fill=0)
rare_syn.idxs <- which(coding.variants$csq == "synonymous"
                              & coding.variants$AF < 0.01)
rare_syn.vids <- unlist(coding.variants$variant_ids[rare_syn.idxs])
rare_syn.ac <- query.ad.matrix(germline.ad, vids=rare_syn.vids,
                                      action="sum", missing.vid.fill=0)

# Get somatic mutation status
Y.vals <- query.ad.matrix(somatic.ad, somatic.vids, elig.controls=elig.controls,
                          action="any", missing.vid.fill=0)
Y.vals <- Y.vals[which(!is.na(Y.vals))]

# Build simple analysis dataframe for subsequent cox modeling
meta.samples <- rownames(meta)
Y.samples <- names(Y.vals)
final.samples <- intersect(meta.samples, Y.samples)
surv.df <- meta[final.samples,
                c("IS_ALIVE", "DAYS_SURVIVED", "SEX",
                  "AGE_AT_DIAGNOSIS", "APPROX_STAGE",
                  colnames(meta)[grep("^cohort.", colnames(meta))])]
surv.df$DEAD <- 1 - surv.df$IS_ALIVE
surv.df$SEX <- relevel(as.factor(surv.df$SEX), ref="FEMALE")
surv.df$KRAS <- Y.vals[final.samples]
surv.df$V10I <- v10i.ac[final.samples]
surv.df$G388R <- g388r.ac[final.samples]
surv.df$rare_nonsyn <- rare_nonsyn.ac[final.samples]
surv.df$rare_syn <- rare_syn.ac[final.samples]
surv.df$ADVANCED_DISEASE <- surv.df$APPROX_STAGE > 2

# Assign sample IDs to each cohort
cohort.sids <- lapply(meta.list, function(l){intersect(rownames(surv.df), l[, 1])})
names(cohort.sids) <- cohort.names

# # Impute missing values on a cohort-specific basis
# for(sids in cohort.sids){
#   surv.df[sids, ] <- impute.missing.values(surv.df[sids, ])
# }

# Get sample IDs per KRAS status
kras.sids <- lapply(1:2, function(ki){
  rownames(surv.df)[which(surv.df$KRAS == ki - 1)]
})
names(kras.sids) <- c("WT", "Mut")

# Fit cohort-specific Cox models
cohort.cox.res <- lapply(1:n.cohorts, function(ci){
  # Further stratify by KRAS mutation status
  c.res <- lapply(1:2, function(ki){
    fit.cox(surv.df[intersect(cohort.sids[[ci]], kras.sids[[ki]]), ],
            other.terms=c("ADVANCED_DISEASE", "V10I", "G388R", "rare_nonsyn", "rare_syn"))
  })
  names(c.res) <- names(kras.sids)
  return(c.res)
})
names(cohort.cox.res) <- cohort.names

# Fit pooled Cox model of all cohorts
pooled.cox.res <- lapply(1:2, function(ki){
  fit.cox(surv.df[kras.sids[[ki]], ],
          other.terms=c("ADVANCED_DISEASE", "V10I", "G388R", "rare_nonsyn", "rare_syn",
                        colnames(surv.df)[grep("^cohort.", colnames(surv.df))]))
})
names(pooled.cox.res) <- names(kras.sids)

# Meta-analyze survival models for each cohort
# TODO: implement this

# Visualize survival vs. KRAS per cohort, split by advanced disease
# TODO: implement this

# Visualize survival vs. FGFR4 G388R + rare nonsyn per cohort among KRAS WT, split by advanced disease
# TODO: implement this

# Visualize survival vs. FGFR4 G388R + rare nonsyn  per cohort among KRAS mutants, split by advanced disease
# TODO: implement this




# Establish baseline KRAS survival difference in PROFILE
kras.wt.ids <- rownames(surv.df)[which(surv.df$KRAS == 0)]
kras.wt <- surv.df[kras.wt.ids, ]
kras.tier1.ids <- rownames(surv.df)[which(surv.df$KRAS == 1)]
kras.tier1 <- surv.df[kras.tier1.ids, ]
km.curve(list("KRAS_WT" = summary(survfit(Surv(kras.wt$DAYS_SURVIVED, 1 - kras.wt$IS_ALIVE) ~ 0)),
              "KRAS_tier1" = summary(survfit(Surv(kras.tier1$DAYS_SURVIVED, 1 - kras.tier1$IS_ALIVE) ~ 0))),
         colors=c("gray50", "red3"), parmar=c(2, 3, 0.25, 8), xlims=c(0, 10*365))

# Establish baseline FGFR4 survival in all samples
fgfr4.g388r.ref.ids <- intersect(rownames(surv.df), names(which(g388r.ac == 0)))
fgfr4.g388r.ref <- surv.df[fgfr4.g388r.ref.ids, ]
fgfr4.g388r.het.ids <- intersect(rownames(surv.df), names(which(g388r.ac == 1)))
fgfr4.g388r.het <- surv.df[fgfr4.g388r.het.ids, ]
fgfr4.g388r.hom.ids <- intersect(rownames(surv.df), names(which(g388r.ac == 2)))
fgfr4.g388r.hom <- surv.df[fgfr4.g388r.hom.ids, ]
fgfr4.rare.ids <- intersect(rownames(surv.df), names(which(rare_nonsyn.ac > 0)))
fgfr4.rare <- meta[fgfr4.rare.ids, ]
pdf("~/scratch/FGFR4.CRC.survival.for_Kevin.pdf", height=3.5, width=5)
km.curve(list("G388R_ref" = summary(survfit(Surv(fgfr4.g388r.ref$DAYS_SURVIVED, 1 - fgfr4.g388r.ref$IS_ALIVE) ~ 0)),
              "G388R_het" = summary(survfit(Surv(fgfr4.g388r.het$DAYS_SURVIVED, 1 - fgfr4.g388r.het$IS_ALIVE) ~ 0)),
              "G388R_hom" = summary(survfit(Surv(fgfr4.g388r.hom$DAYS_SURVIVED, 1 - fgfr4.g388r.hom$IS_ALIVE) ~ 0)),
              "FGFR4_rare" = summary(survfit(Surv(fgfr4.rare$DAYS_SURVIVED, 1 - fgfr4.rare$IS_ALIVE) ~ 0))),
         colors=c(csq.colors, "black"), parmar=c(2, 3, 0.25, 8), xlims=c(0, 10*365))
dev.off()
fit.cox(surv.df, other.terms=c("ADVANCED_DISEASE", "V10I", "G388R", "rare_nonsyn", "rare_syn"))

# Evaluate FGFR4 survival in KRAS WT samples
kras.wt.fgfr4.g388r.ref <- surv.df[intersect(fgfr4.g388r.ref.ids, kras.wt.ids), ]
kras.wt.fgfr4.g388r.het <- surv.df[intersect(fgfr4.g388r.het.ids, kras.wt.ids), ]
kras.wt.fgfr4.g388r.hom <- surv.df[intersect(fgfr4.g388r.hom.ids, kras.wt.ids), ]
kras.wt.fgfr4.rare <- surv.df[intersect(fgfr4.rare.ids, kras.wt.ids), ]
km.curve(list("G388R_ref" = summary(survfit(Surv(kras.wt.fgfr4.g388r.ref$DAYS_SURVIVED, 1 - kras.wt.fgfr4.g388r.ref$IS_ALIVE) ~ 0)),
              "G388R_het" = summary(survfit(Surv(kras.wt.fgfr4.g388r.het$DAYS_SURVIVED, 1 - kras.wt.fgfr4.g388r.het$IS_ALIVE) ~ 0)),
              "G388R_hom" = summary(survfit(Surv(kras.wt.fgfr4.g388r.hom$DAYS_SURVIVED, 1 - kras.wt.fgfr4.g388r.hom$IS_ALIVE) ~ 0)),
              "FGFR4_rare" = summary(survfit(Surv(kras.wt.fgfr4.rare$DAYS_SURVIVED, 1 - kras.wt.fgfr4.rare$IS_ALIVE) ~ 0))),
         colors=c(csq.colors, "black"), parmar=c(2, 3, 0.25, 8), xlims=c(0, 10*365))
#
# Evaluate FGFR4 survival in KRAS tier 1 samples
kras.tier1.fgfr4.g388r.ref <- surv.df[intersect(fgfr4.g388r.ref.ids, kras.tier1.ids), ]
kras.tier1.fgfr4.g388r.het <- surv.df[intersect(fgfr4.g388r.het.ids, kras.tier1.ids), ]
kras.tier1.fgfr4.g388r.hom <- surv.df[intersect(fgfr4.g388r.hom.ids, kras.tier1.ids), ]
kras.tier1.fgfr4.g388r.any <- surv.df[intersect(union(fgfr4.g388r.het.ids, fgfr4.g388r.hom.ids), kras.tier1.ids), ]
kras.tier1.fgfr4.rare <- surv.df[intersect(fgfr4.rare.ids, kras.tier1.ids), ]
km.curve(list("Ref" = summary(survfit(Surv(kras.tier1.fgfr4.g388r.ref$DAYS_SURVIVED, 1 - kras.tier1.fgfr4.g388r.ref$IS_ALIVE) ~ 0)),
              "G388R" = summary(survfit(Surv(kras.tier1.fgfr4.g388r.any$DAYS_SURVIVED, 1 - kras.tier1.fgfr4.g388r.any$IS_ALIVE) ~ 0)),
              "G388R_het" = summary(survfit(Surv(kras.tier1.fgfr4.g388r.het$DAYS_SURVIVED, 1 - kras.tier1.fgfr4.g388r.het$IS_ALIVE) ~ 0)),
              "G388R_hom" = summary(survfit(Surv(kras.tier1.fgfr4.g388r.hom$DAYS_SURVIVED, 1 - kras.tier1.fgfr4.g388r.hom$IS_ALIVE) ~ 0)),
              "Rare" = summary(survfit(Surv(kras.tier1.fgfr4.rare$DAYS_SURVIVED, 1 - kras.tier1.fgfr4.rare$IS_ALIVE) ~ 0))),
         colors=c(csq.colors, "black"), parmar=c(2, 3, 0.25, 8), xlims=c(0, 10*365))

# # Cox models
# cox.df <- merge(surv.df, meta, by="row.names", suffixes=c("", ".redundant"), all=F, sort=F)
# cox.df <- impute.missing.values(cox.df)
# coxph(Surv(DAYS_SURVIVED, 1 - IS_ALIVE) ~ SEX + AGE_AT_DIAGNOSIS + (APPROX_STAGE > 2) + (G388R > 0) + (V10I > 0) + (H3 > 0) + cohort.PROFILE + cohort.HMF + nonsyn + rare_syn, data=cox.df)
# coxph(Surv(DAYS_SURVIVED, 1 - IS_ALIVE) ~ SEX + AGE_AT_DIAGNOSIS + (APPROX_STAGE > 2) + (G388R > 0) + (V10I > 0) + (H3 > 0) + cohort.PROFILE + cohort.HMF + nonsyn + rare_syn, data=cox.df[which(cox.df$Y == 0), ])
# coxph(Surv(DAYS_SURVIVED, 1 - IS_ALIVE) ~ SEX + AGE_AT_DIAGNOSIS + (APPROX_STAGE > 2) + (G388R > 0) + (V10I > 0) + (H3 > 0) + cohort.PROFILE + cohort.HMF + nonsyn + rare_syn, data=cox.df[which(cox.df$Y == 1), ])
