#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Conduct germline-somatic association mega-analyses pooled across multiple cohorts


#########
# Setup #
#########
# Load necessary libraries and constants
require(logistf, quietly=TRUE)
require(argparse, quietly=TRUE)
require(OncoModR, quietly=TRUE)
OncoModR::load.constants("names")


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description=paste("Conduct pooled germline-somatic association",
                                           "tests across two or more cohorts"))
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
parser$add_argument("--somatic-variant-sets", metavar=".tsv", type="character",
                    help="Two-column .tsv of somatic variant sets", required=TRUE)
parser$add_argument("--germline-variant-sets", metavar=".tsv", type="character",
                    help="Two-column .tsv of germline variant sets", required=TRUE)
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="output .tsv file for association statistics")
parser$add_argument("--germline-gq", metavar=".tsv", type="character", action="append",
                    help=paste("Matrix of germline genotype qualities. One ",
                               "per cohort, in order. If provided, GQ will be ",
                               "included as a covariate in association testing. ",
                               "All variants in --somatic-ad must also be present ",
                               "in these files."))
parser$add_argument("--eligible-controls", metavar="path", type="character",
                    help=paste("path to one or more lists of samples eligible to",
                               "be treated as controls for association testing.",
                               "If no file is provided, all samples are eligible."),
                    action="append")
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

# # DEV:
# args <- list("sample_metadata" = c("~/scratch/TCGA.ALL.sample_metadata.tsv.gz",
#                                    "~/scratch/PROFILE.ALL.sample_metadata.tsv.gz",
#                                    "~/scratch/HMF.ALL.sample_metadata.tsv.gz"),
#              "somatic_ad" = c("~/scratch/TCGA.somatic_variants.dosage.tsv.gz",
#                               "~/scratch/PROFILE.somatic_variants.dosage.tsv.gz",
#                               "~/scratch/HMF.somatic_variants.dosage.tsv.gz"),
#              "germline_ad" = c("~/scratch/TCGA.RAS_loci.dosage.sub.tsv.gz",
#                                "~/scratch/PROFILE.RAS_loci.dosage.sub.tsv.gz",
#                                "~/scratch/HMF.RAS_loci.dosage.sub.tsv.gz"),
#              "name" = c("TCGA", "PROFILE", "HMF"),
#              "somatic_variant_sets" = "~/scratch/PDAC.KRAS.somatic_endpoints.tsv",
#              "germline_variant_sets" = "~/scratch/PDAC.KRAS.germline_sets.shard_1",
#              "outfile" = "~/scratch/pooled.assoc.test.tsv",
#              "germline_gq" = c("~/scratch/TCGA.RAS_loci.GQ.sub.tsv.gz",
#                                "~/scratch/PROFILE.RAS_loci.GQ.sub.tsv.gz",
#                                "~/scratch/HMF.RAS_loci.GQ.sub.tsv.gz"),
#              "eligible_controls" = c("~/scratch/TCGA.ALL.eligible_controls.list",
#                                      "~/scratch/PROFILE.ALL.eligible_controls.list",
#                                      "~/scratch/HMF.ALL.eligible_controls.list"),
#              "cancer_type" = "PDAC",
#              "normalize_germline_ad" = FALSE,
#              "multiPop_min_ac" = 10,
#              "multiPop_min_freq" = 0.01)

# # Firth fallback debug DEV:
# Note: not yet updated for GQ adjustment
# args <- list("sample_metadata" = c("~/scratch/TCGA.ALL.sample_metadata.tsv.gz",
#                                    "~/scratch/PROFILE.ALL.sample_metadata.tsv.gz",
#                                    "~/scratch/HMF.ALL.sample_metadata.tsv.gz"),
#              "somatic_ad" = c("~/scratch/TCGA.somatic_variants.dosage.dbg_sub.tsv.gz",
#                               "~/scratch/PROFILE.somatic_variants.dosage.dbg_sub.tsv.gz",
#                               "~/scratch/HMF.somatic_variants.dosage.dbg_sub.tsv.gz"),
#              "germline_ad" = c("~/scratch/TCGA.RAS_loci.dosage.dbg_sub.tsv.gz",
#                                "~/scratch/PROFILE.RAS_loci.dosage.dbg_sub.tsv.gz",
#                                "~/scratch/HMF.RAS_loci.dosage.dbg_sub.tsv.gz"),
#              "name" = c("TCGA", "PROFILE", "HMF"),
#              "somatic_variant_sets" = "~/scratch/dbg.somatic_variants.tsv",
#              "germline_variant_sets" = "~/scratch/dbg.RAS_loci.tsv",
#              "outfile" = "~/scratch/pooled.assoc.test.dbg.tsv",
#              "eligible_controls" = c("~/scratch/TCGA.ALL.eligible_controls.list",
#                                      "~/scratch/PROFILE.ALL.eligible_controls.list",
#                                      "~/scratch/HMF.ALL.eligible_controls.list"),
#              "cancer_type" = "PDAC",
#              "normalize_germline_ad" = FALSE,
#              "multiPop_min_ac" = 10,
#              "multiPop_min_freq" = 0.01)


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
if(length(args$germline_gq) > 0 &
   length(args$germline_gq) != n.cohorts){
  stop("Must either not specify --germline-gq or provide exactly one per cohort. Exiting.")
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
germline.sets <- load.variant.sets(args$germline_variant_sets)
all.germline.vids <- unique(unlist(sapply(germline.sets$variant_ids, parse.variant.set)))
somatic.sets <- load.variant.sets(args$somatic_variant_sets)
all.somatic.vids <- unique(unlist(sapply(somatic.sets$variant_ids, parse.variant.set)))

# Load germline and somatic allele depth matrixes
germline.ad <- lapply(args$germline_ad, load.ad.matrix, sample.subset=samples.w.pheno,
                      variant.subset=all.germline.vids,
                      normalize=args$normalize_germline_ad)
somatic.ad <- lapply(args$somatic_ad, load.ad.matrix, sample.subset=samples.w.pheno,
                     variant.subset=all.somatic.vids)

# Load germline GQ matrix, if provided
if(!is.null(args$germline_gq) & length(args$germline_gq) > 0){
  germline.gq <- lapply(args$germline_gq, load.ad.matrix,
                        sample.subset=samples.w.pheno,
                        variant.subset=all.germline.vids,
                        normalize=FALSE)
}else{
  germline.gq <- NULL
}

# Compute association statistics for each somatic endpoint
res.by.somatic <- apply(somatic.sets, 1, function(somatic.info){
  som.sid <- as.character(somatic.info[1])
  som.vids <- parse.variant.set(somatic.info[2])
  if(!all(sapply(som.vids, function(vids){any(vids %in% unlist(sapply(somatic.ad, rownames)))}))){
    return(NULL)
  }

  # Require at least one each of somatic carriers and reference to proceed
  # Note that we use missing.vid.fill=0 here for somatic mutations because it is
  # safe to assume that (most of the time) effectively all coding KRAS mutations
  # should theoretically be detectable for all samples in all cohorts, and the
  # lack of a reported mutation is not due to lack of data but because no samples
  # carried that mutation (and thus all samples should be reported as AC=0)
  if(is.list(som.vids) & length(som.vids) > 1){
    y.parts <- sapply(som.vids, query.ad.matrix, ad=somatic.ad,
                      elig.controls=elig.controls, action="any",
                      missing.vid.fill=0)
    y.vals <- apply(y.parts, 1, min)
  }else{
    y.vals <- query.ad.matrix(somatic.ad, unlist(som.vids), elig.controls,
                              action="any", missing.vid.fill=0)
  }
  if(length(table(y.vals)) < 2){
    return(NULL)
  }

  # Apply over germline sets
  res.by.germline <- apply(germline.sets, 1, function(germline.info){
    germ.sid <- as.character(germline.info[1])
    germ.vids <- unlist(parse.variant.set(germline.info[2]))
    if(!any(germ.vids %in% unlist(sapply(germline.ad, rownames)))){
      return(NULL)
    }
    x.vals.all <- query.ad.matrix(germline.ad, germ.vids, action="verbose")

    # Prepare GQ data, if --germline-gq is provided
    if(!is.null(germline.gq)){
      gq.vals <- query.gq.matrix(germline.gq, germline.ad, vids=germ.vids)
    }else{
      gq.vals <- NULL
    }

    # Compress x.vals as sum of ACs
    x.vals <- compress.ad.matrix(x.vals.all, action="sum", na.behavior="all")

    # Run germline-somatic association
    res <- tryCatch(germline.somatic.assoc(y.vals, x.vals, meta, gqs=gq.vals,
                                           strict.fallback=F,
                                           custom.covariates=colnames(meta)[grep("^cohort\\.", colnames(meta))],
                                           multiPop.min.ac=args$multiPop_min_ac,
                                           multiPop.min.freq=args$multiPop_min_freq),
                    error=function(e){
                      # If model fails to converge, try again with a simpler model
                      # by dropping lower-rank PCs below PCs 1-3
                      tryCatch(germline.somatic.assoc(y.vals, x.vals,
                                             meta[, grep("^PC[4-9].|^PC[1-9][0-9].", colnames(meta), invert=TRUE)],
                                             gqs=gq.vals, strict.fallback=F,
                                             custom.covariates=colnames(meta)[grep("^cohort\\.", colnames(meta))],
                                             multiPop.min.ac=args$multiPop_min_ac,
                                             multiPop.min.freq=args$multiPop_min_freq,
                                             model.suffix=".simple_PCs"),
                               error=function(e){
                                 # If model with fewer PCs fails to converge, then return null vector
                                 germline.somatic.assoc(y.vals, x.vals,
                                                        meta[, grep("^PC[4-9].|^PC[1-9][0-9].", colnames(meta), invert=TRUE)],
                                                        gqs=gq.vals, strict.fallback=F,
                                                        custom.covariates=colnames(meta)[grep("^cohort\\.", colnames(meta))],
                                                        multiPop.min.ac=args$multiPop_min_ac,
                                                        multiPop.min.freq=args$multiPop_min_freq,
                                                        only.return.freqs=TRUE)
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

# Write cleaned results to --outfile
colnames(res)[1] <- paste("#", colnames(res)[1], sep="")
write.table(res, args$outfile, col.names=T, row.names=F, sep="\t", quote=F)
