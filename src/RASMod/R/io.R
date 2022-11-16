#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Data I/O handlers


#' Load patient metadata
#'
#' Load a standard patient metadata file for a single cohort
#'
#' @param file Path to input patient metadata .tsv
#' @param fill.missing Behavior for handling missing values. See `Details`. \[default: NA\]
#' @param fill.columns Specify in which columns missing values should be filled.
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
  df <- read.table(file, header=T, sep="\t", comment.)
}
