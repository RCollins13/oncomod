#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Basic association testing functions


#' Query allele dosage matrix
#'
#' Subsets and optionally compresses an allele dosage matrix
#'
#' @param ad Allele dosage matrix as imported by [load.ad.matrix]
#' @param vids Variant IDs to query
#' @param action Action to apply to query rows; see `Details`
#'
#' @details Recognized values for `action` include:
#' * `"verbose"` : return the full query matrix \[default\]
#' * `"any"` : return numeric indicator if any of the query rows are non-zero
#' for each sample
#' * `"sum"` : return the sum of allele dosages for all query rows per sample
#'
#' @return numeric vector or data.frame, depending on `action`
#'
#' @export query.ad.matrix
#' @export
query.ad.matrix <- function(ad, vids, action="verbose"){

  # Subset to variants of interest and return directly if optioned
  sub.df <- ad[which(rownames(ad) %in% vids), ]
  if(action == "verbose"){
    return(sub.df)
  }

  # Otherwise, apply various compression strategies prior to returning
  col.all.na <- apply(sub.df, 2, function(vals){all(is.na(vals))})
  if(action == "any"){
    query.res <- as.numeric(apply(sub.df, 2, function(vals){any(as.logical(vals), na.rm=T)}))
  }else if(action == "sum"){
    query.res <- apply(sub.df, 2, sum)
  }
  query.res[col.all.na] <- NA
  names(query.res) <- colnames(sub.df)
  return(query.res)
}


#' Inverse-variance weighted meta-analysis
#'
#' Conduct an inverse-variance weighted meta-analysis of effect sizes
#'
#' @param x vector of observed effect sizes
#' @param se vector of observed standard errors
#'
#' @return named vector of meta-analysis effect size and standard error
#'
#' @export ivw.meta
#' @export
ivw.meta <- function(x, se){
  # Compute weights as
  vars <- se^2
  weights <- 1 / vars

  # Compute weighted mean
  x.bar <- weighted.mean(x, weights, na.rm=T)

  # Compute standard error of sum of weighted observations
  se.bar <- sqrt(1 / sum(weights, na.rm=T))

  return(c("statistic" = x.bar, "std.err" = se.bar))
}


