#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Basic plotting functions to generate entire plots
# See plot.helpers.R for smaller subroutines used to generate certain plot elements


#' Scale-Aware Barplot Cluster
#'
#' Plot a set of stacked barplots scaled proportional to set size
#'
#' @param values List of character vectors of values to plot
#' @param colors Named vector of colors to map to `values`
#' @param group.names (Optional) group names to assign to each list element in `values`
#' @param sep.wex Relative width scalar for whitespace on the X-axis
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
#' @seealso [OncoModR::scaled.swarm], [OncoModR::clean.axis]
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
  group.lefts <- c(0, cumsum(group.widths)[-n.groups]) + (sep.wex * (0:(n.groups-1)))
  group.rights <- group.widths + group.lefts
  group.mids <- (group.lefts + group.rights) / 2
  xlims <- c(-sep.wex, max(group.rights))

  # Prep plot area
  prep.plot.area(xlims, c(1, 0), parmar, xaxs="i", yaxs="i")

  # Add X axis
  x.axis.at <- smart.spacing(group.mids, min.dist=0.15)
  axis(1, at=x.axis.at, line=-1, tick=F, labels=group.names, xpd=T)

  # Add Y axis
  clean.axis(2, at=seq(0, 1, 0.25), labels=rev(paste(seq(0, 100, 25), "%", sep="")),
             cex.axis=5/6, infinite=FALSE, label.line=-0.7, title=NULL)

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


#' Scale-Aware Beeswarm Cluster
#'
#' Plot a set of beeswarm distributions scaled proportional to set size
#'
#' @param values List of numeric vectors of values to plot
#' @param colors Vector of colors for the list elements in `values`
#' @param group.names (Optional) group names to assign to each list element in `values`
#' @param sep.wex Relative width scalar for whitespace on the X-axis
#' between groups \[default: 0.05\]
#' @param pch Value passed to `beeswarm`
#' @param pt.cex Value of `cex` passed to `beeswarm`
#' @param y.title Title of Y-axis
#' @param y.title.line Value of `line` for `y.title`
#' @param y.axis.at Custom Y-axis tick positions, if desired
#' @param y.axis.labels Custom Y-axis tick labels, if desired
#' @param parmar Margin values passed to par()
#'
#' @details If `values` is supplied as a named list, those names will be used as
#' `group.names` unless `group.names` is explicitly specified. Otherwise, `group.names`
#' will be set to the ordinal value of each group in `values`.
#'
#' @seealso [OncoModR::scaled.bars], [OncoModR::clean.axis]
#'
#' @export scaled.swarm
#' @export
scaled.swarm <- function(values, colors, group.names=NULL, sep.wex=0.05,
                         pch=19, pt.cex=0.2, y.title=NULL, y.title.line=0.5,
                         y.axis.at=NULL, y.axis.labels=NULL,
                         parmar=c(1, 2.5, 0.25, 0.25)){
  # Ensure beeswarm library is loaded
  require(beeswarm, quietly=TRUE)

  # Summarize plotting data
  values <- lapply(values, function(v){as.numeric(v[!is.na(v)])})
  if(is.null(group.names)){
    group.names <- names(values)
  }
  drop.groups <- sapply(values, length) == 0
  if(any(drop.groups)){
    drop.group.idx <- which(drop.groups)
    values <- values[-drop.group.idx]
    colors <- colors[-drop.group.idx]
    group.names <- group.names[-drop.group.idx]
  }
  n.groups <- length(values)
  group.size <- sapply(values, length)

  # Get plot dimensions
  group.widths <- group.size / sum(group.size)
  group.lefts <- c(0, group.widths[-n.groups]) + (sep.wex * (0:(n.groups-1)))
  group.rights <- group.widths + group.lefts
  group.mids <- (group.lefts + group.rights) / 2
  xlims <- c(-sep.wex, max(group.rights))
  ylims <- range(unlist(values), na.rm=T)

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="r", yaxs="r")

  # Add X axis
  x.axis.at <- smart.spacing(group.mids, min.dist=0.15)
  axis(1, at=x.axis.at, line=-1, tick=F, labels=group.names, xpd=T)

  # Add Y axis
  clean.axis(2, at=y.axis.at, labels=y.axis.labels,
             cex.axis=5/6, infinite=TRUE, label.line=-0.7,
             title.line=y.title.line, title=y.title)

  # Add boxplots and swarms
  sapply(1:n.groups, function(i){
    if(length(values[[i]] > 0)){
      box.buffer <- min(group.widths[i] / 8, sep.wex / 2)
      metrics <- summary(values[[i]])
      box.x <- (group.lefts + group.mids)[i]/2
      q1q3 <- metrics[c("1st Qu.", "3rd Qu.")]
      segments(x0=box.x, x1=box.x,
               y0=max(metrics["Min."], q1q3[1] - (1.5 * diff(q1q3))),
               y1=min(metrics["Max."], q1q3[2] + (1.5 * diff(q1q3))),
               lwd=3, lend="butt", col=colors[i])
      rect(xleft=group.lefts[i] + box.buffer,
           xright=group.mids[i] - box.buffer,
           ybottom=q1q3[1], ytop=q1q3[2],
           border=NA, bty="n", col=colors[i])
      segments(x0=group.lefts[i],  x1=group.mids[i],
               y0=metrics["Median"], y1=metrics["Median"],
               col="white", lend="butt", lwd=1.5)
      beeswarm(values[[i]], at=group.mids[i], side=1,
               pch=pch, cex=pt.cex, col=colors[i],
               corral="wrap", corralWidth=group.widths[i]/2,
               method="swarm", priority="density", add=T)
    }
  })
}


#' Kaplan-Meier Curves
#'
#' Plot Kaplan-Meier curves for one or more datasets or strata
#'
#' @param surv.models List of one or more [`survival::summary.survfit`] objects
#' @param colors Vector of colors for the list elements in `surv.models`
#' @param group.names (Optional) group names to assign to each list element in `surv.models`
#' @param ci.alpha Transparency value `alpha` for confidence interval shading \[default: 0.15\]
#' @param legend Should a legend be plotted?
#' @param legend.names (Optional) mapping of `values` to labels for legend
#' @param legend.label.spacing Minimum vertical spacing between legend labels \[default: 0.075\]
#' @param title (Optional) Title for plot
#' @param xlims (Optional) two-element vector of start and stop values for X-axis, in days
#' @param parmar Margin values passed to par()
#'
#' @seealso [`survival::Surv`], [`survival::survfit`], [`survival::summary.survfit`]
#'
#' @export km.curve
#' @export
km.curve <- function(surv.models, colors, group.names=NULL, ci.alpha=0.15,
                     legend=TRUE, legend.names=NULL, legend.label.spacing=0.075,
                     title=NULL, xlims=NULL, parmar=c(2, 3, 0.25, 4)){
  # Ensure survival library and OncoModR scale constants are loaded within function scope
  require(survival, quietly=TRUE)
  OncoModR::load.constants("scales", envir=environment())

  # Get plotting values
  if(is.null(group.names)){
    group.names <- names(surv.models)
  }
  n.groups <- length(surv.models)
  if(is.null(legend.names)){
    legend.names <- names(surv.models)
  }
  if(is.null(xlims)){
    xlims <- c(0, max(sapply(surv.models, function(ss){max(ss$time, na.rm=T)})))
  }

  # Prep plot area
  prep.plot.area(xlims, c(0, 1.025), parmar)
  x.ax.step <- max(c(floor(xlims[2] / (365*6)), 1))
  x.ax.years <- seq(0, xlims[2]/365, by=x.ax.step)

  # Add confidence intervals
  # Loop over this twice: first to lay white backgrounds, then add colors
  for(layer in c("white", "colors")){
    sapply(1:n.groups, function(i){
      n.times <- length(surv.models[[i]]$time)
      if(n.times > 1){
        x.bottom <- c(0, OncoModR::stretch.vector(surv.models[[i]]$time, 2)[-2*n.times])
        x.top <- rev(x.bottom)
        y.bottom <- c(1, 1, OncoModR::stretch.vector(surv.models[[i]]$lower, 2)[-c(2*n.times-c(0, 1))])
        y.top <- rev(c(1, 1, OncoModR::stretch.vector(surv.models[[i]]$upper, 2)[-c(2*n.times-c(0, 1))]))
        polygon(x=c(x.bottom, x.top), y=c(y.bottom, y.top), border=NA, bty="n",
                col=if(layer == "white"){"white"}else{adjustcolor(colors[[i]], alpha=ci.alpha)})
      }
    })
  }

  # Add K-M curves
  sapply(1:n.groups, function(i){
    n.times <- length(surv.models[[i]]$time)
    # If summary.survfit returns no data, this is either because
    # there are no patients in this group or nobody died.
    # If the latter, we can plot as a flat line at Y=1 until rmean.endtime (I think?)
    if(surv.models[[i]]$n > 0){
      if(n.times == 0){
        x <- c(0, surv.models[[i]]$rmean.endtime)
        y <- c(1, 1)
      }else{
        x <- c(0, OncoModR::stretch.vector(surv.models[[i]]$time, 2))
        y <- c(1, 1, OncoModR::stretch.vector(surv.models[[i]]$surv, 2))[1:length(x)]
      }
      points(x, y, type="l", col=colors[[i]], lwd=3)
    }
  })

  # Add axes
  clean.axis(1, at=x.ax.years*365, labels=x.ax.years, infinite=TRUE,
             title="Years", label.line=-0.75, title.line=0, tck=-0.0175)
  clean.axis(2, title="Survival Probability", infinite=FALSE, tck=-0.0175)
  mtext(title, side=3, line=0, font=2)

  # Add legend
  if(legend){
    final.y <- sapply(surv.models, function(ss){
      # In the case of no events, ss will have zero rows
      # We can default to Y=1 in this case
      if(length(ss$time) == 0){1}else{
        dist.to.rb <- ss$time - xlims[2]
        if(any(dist.to.rb > 0)){
          closest <- which(dist.to.rb == min(dist.to.rb[which(dist.to.rb >= 0)], na.rm=T) & dist.to.rb >= 0)
          closest <- max(c(1, closest-1))
        }else{
          closest <- which(dist.to.rb == max(dist.to.rb, na.rm=T))
          closest <- max(c(1, closest))
        }
        ss$surv[closest]
      }
    })
    yaxis.legend(legend.names[order(final.y)], x=xlims[2] + (0.05*diff(xlims)),
                 y.positions=final.y[order(final.y)],
                 min.label.spacing=legend.label.spacing,
                 sep.wex=0.05*diff(xlims), colors=colors[order(final.y)])
  }
}

