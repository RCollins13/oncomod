#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Data I/O handlers and basic curation/cleaning steps


#' Load Patient Metadata
#'
#' Load a standard patient metadata file for a single cohort
#'
#' @param file Path to input patient metadata .tsv
#' @param fill.missing Behavior for handling missing values. See `Details`. \[default: NA\]
#' @param fill.columns Specify in which columns missing values should be filled.
#' @param deduplicate Deduplicate based on patient ID / first column \[default: FALSE\]
#' See `Details`. \[default: NULL\]
#'
#' @details See [impute.missing.values] for information on `fill.missing` and `fill.columns`.
#' To disable missing value imputation, provide `fill.missing` = `NA`.
#'
#' @seealso [impute.missing.values]
#'
#' @return data.frame of patient metadata
#'
#' @export load.patient.metadata
#' @export
load.patient.metadata <- function(file, fill.missing=NA, fill.columns=NULL,
                                  deduplicate=FALSE){
  # Read contents of file
  df <- read.table(file, header=T, sep="\t", comment.char="",
                   quote='"', na.strings=c(".", "<NA>"), check.names=F)
  colnames(df)[1] <- gsub("#", "", colnames(df)[1], fixed=T)

  # Deduplicate, if optioned
  if(deduplicate){
    df <- df[!duplicated(df[, 1]), ]
  }

  # Fill missing values, if optioned
  if(!is.na(fill.missing)){
    df <- impute.missing.values(df, fill.missing, fill.columns)
  }

  # Coerce metadata to types as needed
  df$AJCC_STAGE <- as.character(df$AJCC_STAGE)

  return(df)
}


#' Merge patient metadata
#'
#' Merge patient metadata across multiple cohorts
#'
#' @param meta.list List of data frames, one per cohort, as loaded by [load.patient.metadata]
#' @param cancer Subset metadata to a specific cancer type before merging
#' \[default: keep all samples\]
#' @param impute.missing Should missing values be imputed on a cohort specific
#' basis? \[default: TRUE\]
#'
#' @return data.frame
#'
#' @seealso [load.patient.metadata]
#'
#' @export merge.patient.metadata
#' @export
merge.patient.metadata <- function(meta.list, cancer, cohort.names=NULL,
                                   impute.missing=TRUE){
  if(is.null(cohort.names)){
    cohort.names <- paste("cohort", 1:length(meta.list), sep="_")
  }
  n.cohorts <- length(cohort.names)
  for(i in 1:n.cohorts){
    meta.sub <- meta.list[[i]]
    if(!is.null(cancer)){
      meta.sub <- meta.sub[which(meta.sub$CANCER_TYPE == cancer), ]
      cat(paste("Retained", nrow(meta.sub), "samples from cancer type", cancer,
                "in cohort", cohort.names[i],"\n"))
    }
    rownames(meta.sub) <- meta.sub[, 1]
    meta.sub[, 1] <- NULL

    # Impute missing phenotype data as median or mode (depending on variable class)
    if(impute.missing){
      meta.sub <- impute.missing.values(meta.sub, fill.missing="median")
    }

    # Suffix PC covariates with cohort name (to be imputed across cohorts separately)
    pc.idxs <- grep("^PC[0-9]", colnames(meta.sub))
    colnames(meta.sub)[pc.idxs] <- paste(colnames(meta.sub)[pc.idxs], cohort.names[i], sep=".")

    # Add one-hot indicator column for every cohort after the first
    if(i>1){
      meta.sub[, paste("cohort", cohort.names[i], sep=".")] <- 1
    }

    # Update list of per-cohort metadata
    meta.list[[i]] <- meta.sub
  }

  # Combine metadata across cohorts
  meta <- dplyr::bind_rows(meta.list)

  # Fill missing cohort indicators
  cohort.col.idxs <- grep("^cohort\\.", colnames(meta))
  for(cidx in cohort.col.idxs){
    vals <- meta[, cidx]
    na.idxs <- which(is.na(vals))
    if(length(na.idxs) > 0){
      vals[na.idxs] <- 0
    }
    meta[, cidx] <- vals
  }

  return(meta)
}


#' Load variant set membership
#'
#' Load a table linking variant sets to their constituent variant IDs
#'
#' @param file Two-column .tsv mapping variant set IDs \(first column\) to
#' comma-delimited list of constituent variant IDs \(second column\)
#'
#' @return data.frame
#'
#' @export load.variant.sets
#' @export
load.variant.sets <- function(file){
  df <- read.table(file, header=F, sep="\t")
  colnames(df) <- c("set_id", "variant_ids")
  return(df)
}


#' Load an allele dosage matrix
#'
#' Reads an allele dosage matrix into memory with options for subsetting & sorting
#'
#' @param file Tab-delimited matrix of allele dosages for samples \(columns\)
#' by variant IDs \(rows\)
#' @param sort.samples Sort samples axis \[default: TRUE\]
#' @param sample.subset Vector of sample IDs to subset. If not provided, samples
#' will not be subset. \[default: no subsetting\]
#' @param sort.variants Sort variant axis \[default: TRUE\]
#' @param variant.subset Vector of variant IDs to subset. If not provided, variants
#' will not be subset. \[default: no subsetting\]
#' @param drop.ref.records Drop variants from matrix with no observed carriers
#' \[default: TRUE\]
#' @param duplicate.variant.action Action for collapsing rows with identical
#' variant IDs \[default: [max]\]
#' @param normalize Standard normalize each variant across all samples
#' \[default: no normalization\]
#'
#' @return data.frame
#'
#' @export load.ad.matrix
#' @export
load.ad.matrix <- function(file, sort.samples=TRUE, sample.subset=NULL,
                           sort.variants=TRUE, variant.subset=NULL,
                           drop.ref.records=TRUE, duplicate.variant.action=max,
                           normalize=FALSE){
  # Load matrix
  ad <- read.table(file, header=T, sep="\t", comment.char="", check.names=F)

  # Handle duplicate variant IDs
  dup.ids <- unique(ad[, 1][duplicated(ad[, 1])])
  for(dup.id in dup.ids){
    dup.idxs <- which(ad[, 1] == dup.id)
    dedup.values <- apply(ad[dup.idxs, -1], 2, duplicate.variant.action, na.rm=T)
    dedup.values[which(is.infinite(dedup.values))] <- NA
    ad[dup.idxs[1], -1] <- dedup.values
    ad <- ad[-dup.idxs[2:length(dup.idxs)], ]
  }

  # Relocate variant ID to rownames
  rownames(ad) <- ad[, 1]
  ad[, 1] <- NULL

  # Subset & sort samples
  if(!is.null(sample.subset)){
    ad <- ad[, which(colnames(ad) %in% sample.subset)]
  }
  if(sort.samples){
    ad <- ad[, order(colnames(ad))]
  }

  # Drop empty/ref records
  if(drop.ref.records){
    ad <- as.data.frame(ad[apply(ad, 1, sum, na.rm=T) > 0, ])
  }

  # Subset & sort variants
  if(!is.null(variant.subset)){
    ad <- ad[which(rownames(ad) %in% variant.subset), ]
  }
  if(sort.variants){
    ad <- ad[order(rownames(ad)), ]
  }

  # Coerce all values to numeric
  ad[, 1:ncol(ad)] <- apply(ad, 2, as.numeric)

  # Normalize, if optioned
  if(normalize){
    ad[1:nrow(ad), ] <- apply(ad, 1, scale)
  }

  return(ad)
}


#' Load GTF
#'
#' Load features from a GTF file
#'
#' @param gtf.in Path to input GTF file
#' @param regions Optional vector specifying regions to query. See details.
#' @param clean Should GTF be restricted to protein-coding exons and genes only?
#' \[default: FALSE]
#'
#' @details By default, all entries in `gtf.in` will be loaded. If `regions`
#' is specified, [bedr::tabix] will be used to extract query regions of interest.
#' Note that this requires local executibles for `bedtools` and `tabix`. See
#' the `bedr` package documentation for more details.
#'
#' @returns GenomicRanges::GRanges object
#'
#' @export load.gtf
#' @export
load.gtf <- function(gtf.in, regions=NULL, clean=FALSE){
  require(rtracklayer, quietly=TRUE)
  if(is.null(regions)){
    gtf_path <- gtf.in
  }else{
    require(bedr, quietly=TRUE)
    gtf_path <- tempfile()
    write.table(bedr::tabix(regions, gtf.in, check.chr=F), gtf_path,
                sep="\t", col.names=F, row.names=F, quote=F)
  }
  gtf <- rtracklayer::import(gtf_path, format="GTF")
  if(clean){
    gtf[gtf$gene_type == "protein_coding" &
        gtf$type %in% c("exon", "gene")]
  }else{
    gtf
  }
}
