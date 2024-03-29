% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{load.gtf}
\alias{load.gtf}
\title{Load GTF}
\usage{
load.gtf(gtf.in, regions = NULL, clean = FALSE)
}
\arguments{
\item{gtf.in}{Path to input GTF file}

\item{regions}{Optional vector specifying regions to query. See details.}

\item{clean}{Should GTF be restricted to protein-coding exons and genes only?
[default: FALSE]}
}
\value{
GenomicRanges::GRanges object
}
\description{
Load features from a GTF file
}
\details{
By default, all entries in \code{gtf.in} will be loaded. If \code{regions}
is specified, \link[bedr:tabix]{bedr::tabix} will be used to extract query regions of interest.
Note that this requires local executibles for \code{bedtools} and \code{tabix}. See
the \code{bedr} package documentation for more details.
}
