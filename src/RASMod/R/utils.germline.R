#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Utility functions for handling and manipulating germline variant data


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
