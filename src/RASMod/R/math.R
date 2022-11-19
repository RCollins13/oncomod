#!/usr/bin/env R

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


#' Smart Element Spacing
#'
#' Enforce minimum distance between elements in a vector
#'
#' @param ideal.values Numeric values to be spaced
#' @param min.dist Minimum distance between values
#' @param lower.limit Values are not allowed to be placed below this limit
#' \[default: no limits\]
#' @param upper.limit Values are not allowed to be placed above this limit
#' \[default: no limits\]
#'
#' @return Numeric vector
#'
#' @export smart.spacing
#' @export
smart.spacing <- function(ideal.values, min.dist, lower.limit=-Inf, upper.limit=Inf){
  # Curate input values
  val.order <- order(ideal.values)
  n.vals <- length(ideal.values)
  if(n.vals * min.dist > abs(upper.limit - lower.limit)){
    stop("Number of points incompatible with 'min.dist' and specified limits")
  }
  ideal.sorted <- as.numeric(ideal.values[val.order])

  # Enforce lower & upper limits on raw values before spacing
  ideal.sorted[which(ideal.sorted < lower.limit)] <- lower.limit
  ideal.sorted[which(ideal.sorted > upper.limit)] <- upper.limit

  # Iteratively update spacing until best balance is reached
  calc.spacing <- function(ideal.sorted, sig.digits=10){
    round(ideal.sorted[-1] - ideal.sorted[-length(ideal.sorted)], sig.digits)
    }
  spacing <- calc.spacing(ideal.sorted)
  while(any(spacing < min.dist)){
    # Pick closest two points
    # (break ties by taking pair of points with greatest room to be moved)
    spacing.w.limits <- calc.spacing(c(lower.limit, ideal.sorted, upper.limit))
    room.to.move <- sapply(1:length(spacing), function(i){sum(spacing.w.limits[c(i, i+2)])})
    can.move <- which(room.to.move > 0)
    move.priority <- intersect(order(room.to.move, decreasing=TRUE), can.move)
    smallest <- min(spacing)
    smaller.idx <- head(intersect(move.priority, which(spacing == smallest)), 1)
    larger.idx <- smaller.idx + 1

    # Slide points away from each other equally s/t they are exactly min.dist apart
    # Never allow points to be placed beyond limits or jump another point in order
    smaller.max.move <- spacing.w.limits[smaller.idx]
    larger.max.move <- spacing.w.limits[larger.idx+1]
    pad.each <- (min.dist - smallest) / 2
    move.smaller <- min(pad.each, smaller.max.move)
    move.larger <- min(pad.each, larger.max.move)
    ideal.sorted[smaller.idx] <- ideal.sorted[smaller.idx] - move.smaller
    ideal.sorted[larger.idx] <- ideal.sorted[larger.idx] + move.larger

    # Update spacing from ideal.sorted
    spacing <- calc.spacing(ideal.sorted)
  }

  return(ideal.sorted[val.order])
}
