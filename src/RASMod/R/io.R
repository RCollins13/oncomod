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

  return(df)
}


#' Impute missing values
#'
#' Impute missing values into a dataframe
#'
#' @param df Dataframe to impute
#' @param fill.missing Behavior for handling missing values. See `Details`. \[default: "mean"\]
#' @param fill.columns Specify in which columns missing values should be filled.
#'
#' @details Recognized values for `fill.missing` include:
#' * `"mean"` : impute missing values as the mean of all non-missing values
#' * `"median"` : impute missing values as the median of all non-missing values
#' * `"mode"` : impute missing values as the mode of all non-missing values
#'
#' Note that `fill.missing` values of `"mean"` and `"median"` are only applicable to
#' columns with numeric values (e.g., age, BMI, etc.). For non-numeric columns,
#' the behavior of `fill.missing = "mode"` will always be applied.
#'
#' All columns will be subjected to automated missing value imputation unless a
#' non-`NULL` value is supplied to `fill.columns`, in which case only the column
#' names passed to `fill.columns` will be filled.

#' @export impute.missing.values
#' @export
impute.missing.values <- function(df, fill.missing="mean", fill.columns=NULL){
  categorical.action <- RASMod::mode
  if(fill.missing == "mean"){
    numeric.action <- mean
  }else if(fill.missing == "median"){
    numeric.action <- median
  }else if(fill.missing == "mode"){
    numeric.action <- mode
  }

  if(is.null(fill.columns)){
    fill.columns <- colnames(df)[-1]
  }
  for(col in fill.columns){
    na.idxs <- which(is.na(df[, col]))
    if(length(na.idxs) > 0 & length(na.idxs) < nrow(df)){
      if(is.numeric(df[, col])){
        df[na.idxs, col] <- numeric.action(df[-na.idxs, col])
      }else{
        df[na.idxs, col] <- categorical.action(df[-na.idxs, col], break.ties="first")
      }
    }
  }

  return(df)
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
  df$variant_ids <- sapply(df$variant_ids, strsplit, split=",")
  df$variant_ids <- sapply(df$variant_ids, unlist)
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
#'
#' @return data.frame
#'
#' @export load.ad.matrix
#' @export
load.ad.matrix <- function(file, sort.samples=TRUE, sample.subset=NULL,
                           sort.variants=TRUE, variant.subset=NULL,
                           drop.ref.records=TRUE){
  # Load matrix, clean header, and relocate variant ID to rownames
  ad <- read.table(file, header=T, sep="\t", comment.char="", check.names=F)
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

  return(ad)
}

