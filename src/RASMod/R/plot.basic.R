#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Basic plotting functions


#' Prepare Plot Area
#'
#' Prepare a standardized & formatted plot area
#'
#' @param xlims Range of values for X axis
#' @param ylims Range of values for Y axis
#' @param parmar Margin values passed to par()
#' @param xaxs Value of `xaxs` passed to plot()
#' @param yaxs Value of `yaxs` passed to plot()
#'
#' @examples
#' prep.plot.area(xlims=c(0, 5), ylims=(-10, 10), parmar=rep(3, 4));
#'
#' @export prep.plot.area
#' @export
prep.plot.area <- function(xlims, ylims, parmar, xaxs="i", yaxs="i"){
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims, type="n",
       xaxs=xaxs, xlab="", xaxt="n",
       yaxs=yaxs, ylab="", yaxt="n")
}


#' Clean axis
#'
#' Print a clean axis using visually pleasing defaults
#'
#' @param side Value passed to `axis()`. See `?axis` for details.
#' @param at Positions where axis ticks should be plotted \[default: `axTicks(1)`\]
#' @param labels Labels for axis ticks \[default: values of `at`\]
#' @param title Axis title
#' @param tck Value passed to `axis()`. See `?axis` for details. \[default: -0.025\]
#' @param cex.axis Value passed to `axis()`. See `?axis` for details. \[default: 5/6\]
#' @param label.line `line` parameter for axis labels \[default: -0.65\]
#' @param title.line `line` parameter for axis title
#' @param infinite Indicator for the axis to be extended infinitely (without ticks) \[default: FALSE\]
#'
#' @returns NULL
#'
#' @export clean.axis
#' @export
clean.axis <- function(side, at=NULL, labels=NULL, title=NULL, tck=-0.025,
                       cex.axis=5/6, label.line=-0.65, title.line=NULL,
                       infinite=FALSE){
  if(infinite){axis(side, at=c(-10e10, 10e10), tck=0, labels=NA)}
  if(is.null(at)){at <- axTicks(side)}
  if(is.null(labels)){labels <- at}
  if(side %in% c(2, 4)){
    las <- 2
    title.line <- 1
  }else{
    las <- 1
    title.line <- 0.5
  }
  axis(side, at=at, labels=NA, tck=tck)
  sapply(1:length(at), function(i){
    axis(side, at=at[i], labels=labels[i], tick=F, cex.axis=cex.axis,
         las=las, line=label.line)
  })
  axis(side, at=mean(at), tick=F, labels=title, line=title.line)
}


#' Legend-Axis Hybrid
#'
#' Add a legend to the right Y-axis with lines connecting legend labels to
#' specified Y positions
#'
#' @param legend.names Labels to be printed in legend
#' @param x Where the legend will connect to the rest of the plot (in X-axis units)
#' @param y.positions Where should the legend labels be placed (in Y-axis units)
#' @param sep.wex Width expansion term for text relative to `x`
#' @param min.label.spacing Minimum distance between any two labels (in Y-axis units) \[default: 0.1\]
#' @param lower.limit No label will be placed below this value on the Y-axis \[default: no limit\]
#' @param upper.limit No label will be placed above this value on the Y-axis \[default: no limit\]
#' @param colors Line colors connecting labels to plot body \[default: all black\]
#' @param lwd Width of line connecting labels to plot body \[default: 3\]
#'
#' @return NULL
#'
#' @export yaxis.legend
#' @export
yaxis.legend <- function(legend.names, x, y.positions, sep.wex,
                         min.label.spacing=0.1, lower.limit=-Inf,
                         upper.limit=Inf, colors=NULL, lwd=3){
  if(is.null(colors)){
    colors <- "black"
  }
  leg.at <- smart.spacing(y.positions, min.dist=min.label.spacing,
                          lower.limit=lower.limit, upper.limit=upper.limit)
  text(x=x + sep.wex, y=leg.at, labels=legend.names, xpd=T, pos=4)
  segments(x0=x, x1=x + (1.5*sep.wex), y0=y.positions, y1=leg.at,
           lwd=lwd, col=colors, xpd=T, lend="round")
}



#' Scale-Aware Barplot Cluster
#'
#' Plot a set of stacked barplots scaled proportional to set size
#'
#' @param values List of character vectors of values to plot
#' @param colors Named vector of colors to map to `values`
#' @param group.names (Optional) group names to assign to each list element in `values`
#' @param sep.cex Relative width scalar for whitespace on the X-axis
#' between groups \[default: 0.05\]
#' @param title (Optional) title to be printed in top-right corner
#' @param legend Should a legend be plotted?
#' @param legend.names (Optional) mapping of `values` to labels for legend
#' @param parmar Margin values passed to par()
#'
#' @details If `values` is supplied as a named list, those names will be used as
#' `group.names` unless `group.names` is explicitly specified. Otherwise, `group.names`
#' will be set to the ordinal value of each group in `values`.
#'
#' @export scaled.bars
#' @export
scaled.bars <- function(values, colors, group.names=NULL, sep.wex=0.05,
                        title=NULL, legend=TRUE, legend.names=NULL,
                        parmar=c(1, 2, 1.25, 4.5)){
  # Summarize plotting data
  values <- lapply(values, function(v){v[!is.na(v)]})
  elig.vals <- unique(unlist(values))
  if(!all(elig.vals %in% names(colors))){
    stop("Names of 'colors' must cover all values of 'values'")
  }
  colors <- colors[which(names(colors) %in% elig.vals)]
  legend.names <- legend.names[which(names(legend.names) %in% elig.vals)]
  if(is.null(group.names)){
    group.names <- names(values)
  }
  counts <- lapply(values,
                   function(v){sapply(names(colors),
                                      function(k){length(which(v == k))})})
  counts.df <- do.call("cbind", counts)
  n.elig.vals <- length(colors)
  n.groups <- length(values)
  group.size <- apply(counts.df, 2, sum)
  pct.df <- do.call("cbind", lapply(1:n.groups, function(i){as.numeric(counts.df[, i] / group.size[i])}))
  rownames(pct.df) <- names(colors)
  colnames(pct.df) <- group.names
  cumpct.df <- apply(pct.df, 2, cumsum)

  # Get plot dimensions
  group.widths <- group.size / sum(group.size)
  group.lefts <- c(0, group.widths[-n.groups]) + (sep.wex * (0:(n.groups-1)))
  group.rights <- group.widths + group.lefts
  group.mids <- (group.lefts + group.rights) / 2
  xlims <- c(-sep.wex, max(group.rights))

  # Prep plot area
  prep.plot.area(xlims, c(1, 0), parmar, xaxs="i", yaxs="i")

  # Add X axis
  x.axis.at <- smart.spacing(group.mids, min.dist=0.15)
  axis(1, at=x.axis.at, line=-1, tick=F, labels=group.names, xpd=T)

  # Add Y axis
  y.ax.at <- seq(0, 1, 0.25)
  axis(2, at=y.ax.at, tck=-0.025, labels=NA)
  axis(2, at=y.ax.at, tick=F, las=2, cex.axis=5/6, line=-0.7,
       labels=rev(paste(100*y.ax.at, "%", sep="")))

  # Add top scale bar
  scale.k <- floor(log10(sum(group.size)))
  bar.len <- (10^scale.k) / sum(group.size)
  if(bar.len > 0.5){
    scale.k <- scale.k - 1
    bar.len <- (10^scale.k) / sum(group.size)
  }
  if(bar.len < 0.25){
    bar.stretch <- ceiling(0.25 / bar.len)
    bar.len <- bar.len * bar.stretch
  }else{
    bar.stretch <- 1
  }
  segments(x0=mean(xlims) - (0.5 * bar.len),
           x1=mean(xlims) + (0.5 * bar.len),
           y0=-0.03, y1=-0.03, xpd=T, lwd=4, lend="butt")
  text(x=mean(xlims), y=-0.01, pos=3, xpd=T, cex=0.85,
       labels=prettyNum(bar.stretch * (10^scale.k), big.mark=","))

  # Add title
  text(x=xlims[2] - sep.wex, y=-0.05, pos=4, font=2, labels=title, xpd=T)

  # Add legend, if optioned
  if(legend){
    if(is.null(legend.names)){
      legend.names <- names(colors)
    }
    leg.mids <- (c(0, cumpct.df[-n.elig.vals, n.groups]) + cumpct.df[, n.groups])/2
    yaxis.legend(legend.names, xlims[2], leg.mids, sep.wex,
                 min.label.spacing=0.075, lower.limit=0.025, upper.limit=0.975,
                 colors)
  }

  # Add rectangles
  sapply(1:n.groups, function(i){
    rect(xleft=group.lefts[i], xright=group.rights[i],
         ybottom=c(0, cumpct.df[-n.elig.vals, i]), ytop=cumpct.df[, i],
         col=colors, border=NA, bty="n")
    rect(xleft=group.lefts[i], xright=group.rights[i], ybottom=0, ytop=1, col=NA, xpd=T)
  })
}

