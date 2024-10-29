#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Utility functions for handling and manipulating genomic variant data


#' Infer germline position
#'
#' Extract the position of a germline variant based on its variant ID
#'
#' @param vid Germline variant ID
#' @param extract Value to extract, either "position" or "chromosome" \[default: "position"\]
#'
#' @returns Single value (either numeric or character) depending on `extract`
#'
#' @details Note that this function assumes the variant adheres to the
#' ID format of chrom_pos_ref_alt used throughout this study
#'
#' @export infer.germline.position.from.id
#' @export
infer.germline.position.from.id <- function(vid, extract="position"){
  v <- unlist(strsplit(unlist(strsplit(vid, split=";"))[1], split="_"))
  if(extract == "position"){
    as.numeric(v[length(v) - 2])
  }else if(extract == "chromosome"){
    as.character(v[length(v) - 3])
  }
}


#' Simplify allele dosage variant IDs
#'
#' Simplify variant IDs for a single allele dosage matrix
#'
#' @param ad.df Allele dosage matrix as read by [RLCtools::query.ad.matrix()]
#'
#' @returns Allele dosage matrix with simplified variant IDs.
#'
#' @details Note that this function assumes the variant IDs adhere to the
#' ID format of chrom_pos_ref_alt used throughout this study
#'
#' @export simplify.vids
#' @export
simplify.vids <- function(ad.df){
  rownames(ad.df) <- sapply(rownames(ad.df), function(vid){
    parts <- unlist(strsplit(vid, split="_"))
    paste(parts[(length(parts)-3):length(parts)], collapse="_")
  })
  return(ad.df)
}

