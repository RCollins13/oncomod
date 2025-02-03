#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Conduct germline-somatic association testing for a single cancer type & cohort


#########
# Setup #
#########
# Load necessary libraries and constants
require(logistf, quietly=TRUE)
require(argparse, quietly=TRUE)
require(OncoModR, quietly=TRUE)
OncoModR::load.constants("names")


##################
# Data Functions #
##################
# Subset a variant set dataframe based on variants present in an AD matrix
subset.variant.sets <- function(variant.sets, ad.matrix){
  keep.idxs <- which(sapply(variant.sets[, 2], function(vid.str){
    any(unlist(parse.variant.set(vid.str)) %in% rownames(ad.matrix))
  }))
  variant.sets[keep.idxs, ]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Conduct single-cohort germline-somatic",
                                           "association tests for a single cancer type"))
parser$add_argument("--sample-metadata", metavar=".tsv", type="character",
                    help="sample metadata .tsv", required=TRUE)
parser$add_argument("--somatic-ad", metavar=".tsv", type="character",
                    help="Somatic allele dosage matrix", required=TRUE)
parser$add_argument("--germline-ad", metavar=".tsv", type="character",
                    help="Germline allele dosage matrix", required=TRUE)
parser$add_argument("--somatic-variant-sets", metavar=".tsv", type="character",
                    help="Two-column .tsv of somatic variant sets", required=TRUE)
parser$add_argument("--germline-variant-sets", metavar=".tsv", type="character",
                    help="Two-column .tsv of germline variant sets", required=TRUE)
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="output .tsv file for association statistics")
parser$add_argument("--germline-gq", metavar=".tsv", type="character",
                    help=paste("Matrix of germline genotype qualities. If ",
                               "provided, GQ will be included as a covariate ",
                               "in association testing. All variants in ",
                               "--somatic-ad must also be present in this file."))
parser$add_argument("--eligible-controls", metavar="path", type="character",
                    help=paste("path to list of samples eligible to be treated as",
                               "controls for association testing. If no file is",
                               "provided, all samples are eligible."))
parser$add_argument("--cancer-type", metavar="character",
                    help=paste("Subset to samples from this cancer type",
                               "[default: use all samples]"))
parser$add_argument("--normalize-germline-ad", action="store_true", default=FALSE,
                    help=paste("Standard normalize each entry in --germline-ad",
                               "before analysis. [default: FALSE]"))
parser$add_argument("--multiPop-min-ac", metavar="integer", default=10, type="integer",
                    help=paste("Restrict tests involving germline or somatic ",
                               "counts below this threshold to European-only ",
                               "[default: 10]"))
parser$add_argument("--multiPop-min-freq", metavar="float", default=0.01, type="double",
                    help=paste("Restrict tests involving germline or somatic ",
                               "frequencies below this threshold to European-only ",
                               "[default: 0.01]"))
args <- parser$parse_args()

# # DEV - TCGA
# args <- list("sample_metadata" = "~/scratch/TCGA.ALL.sample_metadata.tsv.gz",
#              "cancer_type" = "PDAC",
#              "germline_ad" = "~/scratch/TCGA.RAS_loci.dosage.tsv.gz",
#              "germline_variant_sets" = "~/scratch/PDAC.KRAS.germline_sets.shard_1",
#              "germline_gq" = "~/scratch/TCGA.RAS_loci.GQ.tsv.gz",
#              "outfile" = "~/scratch/TCGA.PDAC.KRAS.sumstats.tsv",
#              "somatic_ad" = "~/scratch/TCGA.somatic_variants.dosage.tsv.gz",
#              "somatic_variant_sets" = "~/scratch/PDAC.KRAS.somatic_endpoints.tsv",
#              "eligible_controls" = "~/scratch/TCGA.ALL.eligible_controls.list",
#              "normalize_germline_ad" = FALSE,
#              "multiPop_min_ac" = 10,
#              "multiPop_min_freq" = 0.01)

# # DEV - PROFILE
# # Note: not yet updated for GQ
# args <- list("sample_metadata" = "~/scratch/PROFILE.ALL.sample_metadata.tsv.gz",
#              "cancer_type" = "PDAC",
#              "germline_ad" = "~/scratch/PROFILE.RAS_loci.dosage.tsv.gz",
#              "germline_variant_sets" = "~/scratch/PDAC.NRAS.germline_sets.shard_2",
#              "outfile" = "~/scratch/PROFILE.PDAC.NRAS.sumstats.tsv",
#              "somatic_ad" = "~/scratch/PROFILE.somatic_variants.dosage.tsv.gz",
#              "somatic_variant_sets" = "~/scratch/PDAC.NRAS.somatic_endpoints.tsv",
#              "eligible_controls" = "~/scratch/PROFILE.ALL.eligible_controls.list",
#              "normalize_germline_ad" = FALSE,
#              "multiPop_min_ac" = 10,
#              "multiPop_min_freq" = 0.01)

# # DEV - HMF
# args <- list("sample_metadata" = "~/scratch/HMF.ALL.sample_metadata.tsv.gz",
#              "cancer_type" = "PDAC",
#              "germline_ad" = "~/scratch/HMF.RAS_loci.dosage.dbg.tsv.gz",
#              "germline_variant_sets" = "~/scratch/germline_sets.dbg.tsv",
#              "outfile" = "~/scratch/HMF.PDAC.KRAS.sumstats.tsv",
#              "somatic_ad" = "~/scratch/HMF.somatic_variants.dosage.tsv.gz",
#              "somatic_variant_sets" = "~/scratch/PDAC.KRAS.somatic_endpoints.tsv",
#              "germline_gq" = "~/scratch/HMF.RAS_loci.GQ.dbg.tsv.gz",
#              "eligible_controls" = "~/scratch/HMF.ALL.eligible_controls.list",
#              "normalize_germline_ad" = FALSE,
#              "multiPop_min_ac" = 10,
#              "multiPop_min_freq" = 0.01)

# # DEV - PROFILE PRS
# args <- list("sample_metadata" = "~/scratch/PROFILE.ALL.sample_metadata.tsv.gz",
#              "cancer_type" = "CRAD",
#              "germline_ad" = "~/scratch/PROFILE.PRS.tsv.gz",
#              "germline_variant_sets" = "~/scratch/CRAD.germline_PRS.tsv",
#              "outfile" = "~/scratch/PROFILE.CRAD.KRAS.PRS.sumstats.tsv",
#              "somatic_ad" = "~/scratch/PROFILE.somatic_variants.dosage.sub.tsv.gz",
#              "somatic_variant_sets" = "~/scratch/CRAD.KRAS.somatic_endpoints.tsv",
#              "eligible_controls" = NULL,
#              "normalize_germline_ad" = TRUE,
#              "multiPop_min_ac" = 10,
#              "multiPop_min_freq" = 0.01)

# Load patient metadata and subset to cancer type of interest (if optioned)
meta <- load.patient.metadata(args$sample_metadata, deduplicate=TRUE)
if(!is.null(args$cancer_type)){
  meta <- meta[which(meta$CANCER_TYPE == args$cancer_type), ]
  cat(paste("Retained", nrow(meta), "samples from cancer type", args$cancer_type, "\n"))
}
rownames(meta) <- meta[, 1]
meta[, 1] <- NULL
samples.w.pheno <- rownames(meta)
if(is.null(args$eligible_controls)){
  elig.controls <- samples.w.pheno
}else{
  all.elig.controls <- unique(read.table(args$eligible_controls, header=F)[, 1])
  elig.controls <- intersect(samples.w.pheno, all.elig.controls)
}

# Impute missing phenotype data as median or mode (depending on variable class)
meta <- impute.missing.values(meta, fill.missing="median")

# Load germline and somatic variant lists
germline.sets <- load.variant.sets(args$germline_variant_sets)
all.germline.vids <- unique(unlist(sapply(germline.sets$variant_ids, parse.variant.set)))
somatic.sets <- load.variant.sets(args$somatic_variant_sets)
all.somatic.vids <- unique(unlist(sapply(somatic.sets$variant_ids, parse.variant.set)))

# Load germline and somatic allele depth matrixes
germline.ad <- load.ad.matrix(args$germline_ad, sample.subset=samples.w.pheno,
                              variant.subset=all.germline.vids,
                              normalize=args$normalize_germline_ad)
somatic.ad <- load.ad.matrix(args$somatic_ad, sample.subset=samples.w.pheno,
                             variant.subset=all.somatic.vids)

# Subset variant sets to ADs to avoid iterating over nonexistant germline-somatic pairs
germline.sets <- subset.variant.sets(germline.sets, germline.ad)
somatic.sets <- subset.variant.sets(somatic.sets, somatic.ad)

# Load germline GQ matrix, if provided
if(!is.null(args$germline_gq)){
  germline.gq <- load.ad.matrix(args$germline_gq, sample.subset=samples.w.pheno,
                                variant.subset=all.germline.vids,
                                normalize=FALSE)
}else{
  germline.gq <- NULL
}

# Compute association statistics for each somatic endpoint
res.by.somatic <- apply(somatic.sets, 1, function(somatic.info){
  som.sid <- as.character(somatic.info[1])
  som.vids <- parse.variant.set(somatic.info[2])
  if(!all(sapply(som.vids, function(vids){any(vids %in% rownames(somatic.ad))}))){
    return(NULL)
  }

  # Require at least one each of somatic carriers and reference to proceed
  if(is.list(som.vids) & length(som.vids) > 1){
    y.parts <- sapply(som.vids, query.ad.matrix, ad=somatic.ad,
                      elig.controls=elig.controls, action="any")
    y.vals <- apply(y.parts, 1, min)
  }else{
    y.vals <- query.ad.matrix(somatic.ad, unlist(som.vids), elig.controls, action="any")
  }
  if(length(table(y.vals)) < 2){
    return(NULL)
  }

  # Apply over germline sets
  res.by.germline <- apply(germline.sets, 1, function(germline.info){
    germ.sid <- as.character(germline.info[1])
    germ.vids <- unlist(parse.variant.set(germline.info[2]))
    if(!any(germ.vids %in% rownames(germline.ad))){
      return(NULL)
    }
    x.vals.all <- query.ad.matrix(germline.ad, germ.vids, action="verbose")

    # Prepare GQ data, if --germline-gq is provided
    if(!is.null(germline.gq)){
      gq.vals <- query.gq.matrix(germline.gq, x.vals.all)
    }else{
      gq.vals <- NULL
    }

    # Compress x.vals as sum of ACs
    x.vals <- compress.ad.matrix(x.vals.all, action="sum", na.behavior="all")

    # Prep function aliases with progressively simpler covariates
    m1 <- function(y, x, m){germline.somatic.assoc(y, x, m, gqs=gq.vals,
                                                   multiPop.min.ac=args$multiPop_min_ac,
                                                   multiPop.min.freq=args$multiPop_min_freq)}
    # If model fails to converge, try again with a simpler model
    # by dropping lower-rank PCs below PCs 1-3
    m2 <- function(y, x, m){germline.somatic.assoc(y, x,
                                                   m[, grep("^PC[4-9]|^PC[1-9][0-9]", colnames(meta), invert=TRUE)],
                                                   gqs=gq.vals,
                                                   multiPop.min.ac=args$multiPop_min_ac,
                                                   multiPop.min.freq=args$multiPop_min_freq,
                                                   model.suffix=".simple_PCs")}
    # If model with fewer PCs fails to converge, try
    # one last time by also dropping GQ adjustment
    m3 <- function(y, x, m){germline.somatic.assoc(y, x,
                                                   m[, grep("^PC[4-9]|^PC[1-9][0-9]", colnames(meta), invert=TRUE)],
                                                   multiPop.min.ac=args$multiPop_min_ac,
                                                   multiPop.min.freq=args$multiPop_min_freq,
                                                   model.suffix=".simple_PCs.no_GQ")}
    # If model with fewer PCs and no GQ adjustment
    # fails to converge, then return null vector
    m4 <- function(y, x, m){germline.somatic.assoc(y.vals, x.vals, m,
                                                   strict.fallback=F,
                                                   custom.covariates=colnames(meta)[grep("^cohort\\.", colnames(meta))],
                                                   multiPop.min.ac=args$multiPop_min_ac,
                                                   multiPop.min.freq=args$multiPop_min_freq,
                                                   only.return.freqs=TRUE)}

    # Run germline-somatic association
    res <- tryCatch(m1(y.vals, x.vals, meta),
                    error=function(e){
                      tryCatch(m2(y.vals, x.vals, meta),
                               error=function(e){
                                 tryCatch(m3(y.vals, x.vals, meta),
                                          error=function(e){
                                            m4(y.vals, x.vals, meta)
                                          })
                               })
                    })
    if(!is.null(res)){
      return(c("germline"=germ.sid, res))
    }else{
      return(NULL)
    }
  })
  res.type <- typeof(res.by.germline)
  if(res.type == "list"){
    res.by.germline <- as.data.frame(do.call("rbind", res.by.germline))
  }else if(res.type == "character"){
    res.by.germline <- as.data.frame(t(res.by.germline))
  }
  if(!is.null(res.by.germline)){
    cbind("somatic"=som.sid, res.by.germline)
  }else{
    return(NULL)
  }
})
res <- as.data.frame(do.call("rbind", res.by.somatic))

# Clean up and write results to output file
colnames(res)[1] <- paste("#", colnames(res)[1], sep="")
write.table(res, args$outfile, sep="\t", quote=F, row.names=F, col.names=T)
