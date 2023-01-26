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
# Data functions #
##################
# Run association test for germline & somatic interactions
germline.somatic.assoc <- function(y.vals, x.vals, samples, meta, firth=TRUE){
  # Construct test df from y, x, and meta
  test.df <- meta[samples, c("AGE_AT_DIAGNOSIS", "SEX", paste("PC", 1:10, sep=""),
                             "TUMOR_PURITY")]
  test.df$SEX <- as.numeric(test.df$SEX == "MALE")
  test.df <- cbind(data.frame("Y" = y.vals, "X" = x.vals, row.names=names(y.vals)),
                   test.df)

  # Ensure enough data are still present for association testing
  n.samples <- length(samples)
  somatic.ac <- sum(y.vals, na.rm=T)
  germline.ac <- sum(x.vals, na.rm=T)
  yes_somatic.germline.ac <- sum(x.vals[which(y.vals > 0)], na.rm=T)
  no_somatic.germline.ac <- sum(x.vals[which(y.vals == 0)], na.rm=T)
  if(any(c(n.samples, somatic.ac, germline.ac) == 0)){
    return(c("samples"=n.samples,
             "somatic_AC"=somatic.ac,
             "yes_somatic.germline_AC"=yes_somatic.germline.ac,
             "no_somatic.germline_AC"=no_somatic.germline.ac,
             "beta"=NA,
             "beta_SE"=NA,
             "p"=NA))
  }

  # Drop covariates with no informative observations
  all.na <- which(apply(test.df, 2, function(vals){all(is.na(vals))}))
  if(length(all.na > 0)){
    test.df <- test.df[, -all.na]
  }

  # Standard normalize all covariates
  test.df[, -c(1:2)] <- apply(test.df[, -c(1:2)], 2, scale)

  # Run association model
  if(firth){
    fit <- logistf(Y ~ . + (AGE_AT_DIAGNOSIS * SEX), data=test.df,
                   control=logistf.control(maxit=100))
  }else{
    fit <- glm(Y ~ . + (AGE_AT_DIAGNOSIS * SEX), data=test.df, family="binomial")
  }


  # Extract association stats for germline variants
  if(is.na(fit$coefficients["X"])){
    assoc.res <- rep(NA, 4)
  }else{
    if(firth){
      assoc.res <- as.numeric(c(fit$coefficients["X"],
                              sqrt(diag(vcov(fit)))["X"],
                              qchisq(1-fit$prob, df=1)["X"],
                              fit$prob["X"]))
      stat.name <- "chisq"
    }else{
      assoc.res <- as.numeric(summary(fit)$coefficients["X", ])
      stat.name <- "z"
    }
  }
  # Return summary vector
  res <- c("samples"=n.samples,
    "somatic_AC"=somatic.ac,
    "yes_somatic.germline_AC"=yes_somatic.germline.ac,
    "no_somatic.germline_AC"=no_somatic.germline.ac,
    "beta"=assoc.res[1],
    "beta_SE"=assoc.res[2],
    "statistic"=assoc.res[3],
    "p"=assoc.res[4])
  names(res)[which(names(res) == "statistic")] <- stat.name
  return(res)
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
parser$add_argument("--multiPop-min-ac", metavar="integer", default=10, type="integer",
                    help=paste("Restrict tests involving germline or somatic ",
                               "counts below this threshold to European-only ",
                               "[default: 10]"))
parser$add_argument("--multiPop-min-freq", metavar="float", default=0.01, type="double",
                    help=paste("Restrict tests involving germline or somatic ",
                               "frequencies below this threshold to European-only ",
                               "[default: 0.01]"))
args <- parser$parse_args()
#
# # DEV - TCGA
# args <- list("sample_metadata" = "~/scratch/TCGA.ALL.sample_metadata.tsv.gz",
#              "cancer_type" = "PDAC",
#              "germline_ad" = "~/scratch/TCGA.RAS_loci.dosage.tsv.gz",
#              "germline_variant_sets" = "~/scratch/PDAC.KRAS.germline_sets.tsv",
#              "outfile" = "~/scratch/TCGA.PDAC.KRAS.sumstats.tsv",
#              "somatic_ad" = "~/scratch/TCGA.somatic_variants.dosage.tsv.gz",
#              "somatic_variant_sets" = "~/scratch/PDAC.KRAS.somatic_endpoints.tsv",
#              "multiPop_min_ac" = 10,
#              "multiPop_min_freq" = 0.01)

# # DEV - PROFILE
# args <- list("sample_metadata" = "~/scratch/PROFILE.ALL.sample_metadata.tsv.gz",
#              "cancer_type" = "PDAC",
#              "germline_ad" = "~/scratch/PROFILE.RAS_loci.dosage.tsv.gz",
#              "germline_variant_sets" = "~/scratch/PDAC.HRAS.germline_sets.tsv",
#              "outfile" = "~/scratch/PROFILE.PDAC.HRAS.sumstats.tsv",
#              "somatic_ad" = "~/scratch/PROFILE.somatic_variants.dosage.tsv.gz",
#              "somatic_variant_sets" = "~/scratch/PDAC.HRAS.somatic_endpoints.tsv",
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
                              variant.subset=germline.vids)
somatic.ad <- load.ad.matrix(args$somatic_ad, sample.subset=samples.w.pheno,
                             variant.subset=somatic.vids)

# Compute association statistics for each somatic endpoint
res.by.somatic <- apply(somatic.sets, 1, function(somatic.info){
  som.sid <- as.character(somatic.info[1])
  som.vids <- as.vector(unlist(somatic.info[2]))
  if(!any(som.vids %in% rownames(somatic.ad))){
    return(NULL)
  }
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

    # Enforce strict intersection of non-NA samples between germline & somatic
    som.samples <- names(y.vals)[which(!is.na(y.vals))]
    germ.samples <- names(x.vals)[which(!is.na(x.vals))]
    samples <- intersect(som.samples, germ.samples)
    y.vals <- y.vals[samples]
    x.vals <- x.vals[samples]

    # If germline or somatic variants are <1% freq / <10 counts, restrict to
    # European samples to protect against pop strat
    eur.only <- FALSE
    if(sum(x.vals) < args$multiPop_min_ac
       | sum(x.vals) / length(samples) < args$multiPop_min_freq
       | sum(y.vals) < args$multiPop_min_ac
       | sum(y.vals) / length(samples) < args$multiPop_min_freq){
      eur.samples <- rownames(meta[which(meta$POPULATION == "EUR"), ])
      samples <- intersect(samples, eur.samples)
      x.vals <- x.vals[samples]
      y.vals <- y.vals[samples]
      eur.only <- TRUE
    }

    # Ensure non-zero counts for X and Y
    if(length(samples) == 0
       | length(table(x.vals)) < 2
       | length(table(y.vals)) < 2){
      return(NULL)
    }

    # If all checks pass above, run association test
    c("germline"=germ.sid, germline.somatic.assoc(y.vals, x.vals, samples, meta),
      "EUR_only"=eur.only)
  })
  if(typeof(res.by.germline) == "list"){
    res.by.germline <- as.data.frame(do.call("rbind", res.by.germline))
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
