#!/usr/bin/env bash

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
#' @param
#' See `Details`. \[default: NULL\]
#'
#' @details Recognized values for `fill.missing` include:
#' * `NA` : do not fill any missing values
#' * `"mean"` : impute missing values as the mean of all non-missing values
#' * `"median"` : impute missing values as the median of all non-missing values
#' * `"mode"` : impute missing values as the mode of all non-missing values
#'
#' Note that `fill.missing` values of `"mean"` and `"median"` are only applicable to
#' columns with numeric values (e.g., age, BMI, etc.). For non-numeric columns,
#' the behavior of `fill.missing = "mode"` will always be applied unless `fill.missing`
#' is set to `NA` (the default\).
#'
#' All columns will be subjected to automated missing value imputation unless a
#' non-`NULL` value is supplied to `fill.columns`, in which case only the column
#' names passed to `fill.columns`will be filled.
#'
#' @return data.frame of patient metadata
#'
#' @export load.patient.metadata
#' @export
load.patient.metadata <- function(file, fill.missing=NA, missing.columns=NULL){
  # Read contents of file
  df <- read.table(file, header=T, sep="\t", comment.char="",
                   na.strings=c(".", "<NA>"), check.names=F)
  colnames(df)[1] <- gsub("#", "", colnames(df)[1], fixed=T)

  # Fill missing values, if optioned
  if(!is.na(fill.missing)){
    categorical.action <- RASMod::mode
    if(fill.missing == "mean"){
      numeric.action <- mean
    }else if(fill.missing == "median"){
      numeric.action <- median
    }else if(fill.missing == "mode"){
      numeric.action <- mode
    }

    if(is.null(missing.columns)){
      missing.columns <- colnames(df)[-1]
    }
    for(col in missing.columns){
      na.idxs <- which(is.na(df[, col]))
      if(length(na.idxs) > 0 & length(na.idxs) < nrow(df)){
        if(is.numeric(df[, col])){
          col[na.idxs] <- numeric.action(df[-na.idxs, col])
        }else{
          col[na.idxs] <- categorical.action(df[-na.idxs, col], break.ties="first")
        }
      }
    }
  }

  return(df)
}
