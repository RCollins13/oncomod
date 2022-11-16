#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Mathematical helper functions


#' Arithmetic Mode
#'
#' Compute the arithmetic mode of a vector
#'
#' @param values Vector of values
#' @param break.ties Specify how ties should be broken. See `Details`. \[default: "all"\]
#' @param seed Random seed to break ties, if needed. \[default: 2022\]
#'
#' @details Recognized values for `break.ties` include:
#' * `all` : return all tied modes
#' * `first` : return the first mode value encountered
#' * `last` : return the last mode value encountered
#' * `random` : return a randomly selected mode value
#'
#' @return mode(s) of `values`
#'
#' @export mode
#' @export
mode <- function(values, break.ties="all", seed=2022){
  vt <- sort(table(values), decreasing=TRUE)
  vmax <- max(vt, na.rm=TRUE)
  modes <- names(vt)[which(vt == vmax)]
  if(break.ties == "all"){
    return(modes)
  }else if(break.ties == "first"){
    return(head(modes, 1))
  }else if(break.ties == "last"){
    return(tail(modes, 1))
  }else if(break.ties == "random"){
    set.seed(seed)
    return(sample(modes, 1))
  }
}
