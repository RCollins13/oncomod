#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Miscellaneous helper functions


#' Stretch Vector
#'
#' Expand the length of a vector by repeating (or "stuttering") values
#'
#' @param values Vector of values to be stretched
#' @param k Number of times to duplicate each element in `values`
#'
#' @examples
#' stretch.vector(values=c(1, 2, 3), k=4)
#'
#' @export stretch.vector
#' @export
stretch.vector <- function(values, k){
  as.vector(sapply(values, rep, times=k))
}


#' Parse Variant Set Members
#'
#' Parse a string of variant IDs that comprise a variant set
#'
#' @param vid.str Character object to be parsed
#'
#' @details Variant IDs should be separated by commas.
#'
#' Pipe delimiters will be interpreted as multiple independent sets of variant IDs
#' \(as is the case for co-mutation analyses, etc.\)
#'
#' @return List of character vectors
#'
#' @export parse.variant.set
#' @export
parse.variant.set <- function(vid.str){
  lapply(strsplit(as.character(vid.str), split="|", fixed=T),
                 strsplit, split=",", fixed=T)[[1]]
}

