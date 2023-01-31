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
require(RASMod, quietly=TRUE)
require(logistf, quietly=TRUE)
require(argparse, quietly=TRUE)
RASMod::load.constants("names")


##################
# Data Functions #
##################
# Subset a variant set dataframe based on variants present in an AD matrix
subset.variant.sets <- function(variant.sets, ad.matrix){
  variant.sets[which(unlist(sapply(variant.sets[, 2], function(vids){
    any(vids %in% rownames(ad.matrix))
  }))), ]
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
#              "germline_ad" = "~/scratch/TCGA.RAS_loci.dosage.sub.tsv.gz",
#              "germline_variant_sets" = "~/scratch/PDAC.KRAS.germline_sets.tsv",
#              "outfile" = "~/scratch/TCGA.PDAC.KRAS.sumstats.tsv",
#              "somatic_ad" = "~/scratch/TCGA.somatic_variants.dosage.sub.tsv.gz",
#              "somatic_variant_sets" = "~/scratch/PDAC.KRAS.somatic_endpoints.tsv",
#              "normalize_germline_ad" = FALSE,
#              "multiPop_min_ac" = 10,
#              "multiPop_min_freq" = 0.01)

# # DEV - PROFILE
# args <- list("sample_metadata" = "~/scratch/PROFILE.ALL.sample_metadata.tsv.gz",
#              "cancer_type" = "SKCM",
#              "germline_ad" = "~/scratch/PROFILE.RAS_loci.dosage.tsv.gz",
#              "germline_variant_sets" = "~/scratch/SKCM.KRAS.germline_sets.shard_39",
#              "outfile" = "~/scratch/PROFILE.SKCM.KRAS.sumstats.tsv",
#              "somatic_ad" = "~/scratch/PROFILE.somatic_variants.dosage.tsv.gz",
#              "somatic_variant_sets" = "~/scratch/SKCM.KRAS.somatic_endpoints.tsv",
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

# Impute missing phenotype data as median or mode (depending on variable class)
meta <- impute.missing.values(meta, fill.missing="median")

# Load germline and somatic variant lists
germline.sets <- load.variant.sets(args$germline_variant_sets)
germline.vids <- unique(unlist(germline.sets$variant_ids))
somatic.sets <- load.variant.sets(args$somatic_variant_sets)
somatic.vids <- unique(unlist(somatic.sets$variant_ids))

# Load germline and somatic allele depth matrixes
germline.ad <- load.ad.matrix(args$germline_ad, sample.subset=samples.w.pheno,
                              variant.subset=germline.vids,
                              normalize=args$normalize_germline_ad)
somatic.ad <- load.ad.matrix(args$somatic_ad, sample.subset=samples.w.pheno,
                             variant.subset=somatic.vids)

# Subset variant sets to ADs to avoid iterating over nonexistant germline-somatic pairs
germline.sets <- subset.variant.sets(germline.sets, germline.ad)
somatic.sets <- subset.variant.sets(somatic.sets, somatic.ad)

# Compute association statistics for each somatic endpoint
res.by.somatic <- apply(somatic.sets, 1, function(somatic.info){
  som.sid <- as.character(somatic.info[1])
  som.vids <- as.vector(unlist(somatic.info[2]))
  if(!any(som.vids %in% rownames(somatic.ad))){
    return(NULL)
  }

  # Require at least one each of somatic carriers and reference to proceed
  y.vals <- query.ad.matrix(somatic.ad, som.vids, action="any")
  if(length(table(y.vals)) < 2){
    return(NULL)
  }

  # Apply over germline sets
  res.by.germline <- apply(germline.sets, 1, function(germline.info){
    germ.sid <- as.character(germline.info[1])
    germ.vids <- as.vector(unlist(germline.info[2]))
    if(!any(germ.vids %in% rownames(germline.ad))){
      return(NULL)
    }
    x.vals <- query.ad.matrix(germline.ad, germ.vids, action="sum")

    # Run germline-somatic association
    res <- germline.somatic.assoc(y.vals, x.vals, samples, meta,
                                           multiPop.min.ac=args$multiPop_min_ac,
                                           multiPop.min.freq=args$multiPop_min_freq)
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
