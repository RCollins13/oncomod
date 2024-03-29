% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assoc.meta.R
\name{load.assoc.stats.single}
\alias{load.assoc.stats.single}
\title{Load association stats}
\usage{
load.assoc.stats.single(file, suffix = NULL, suffix.skip = 1:2)
}
\arguments{
\item{file}{Path to input stats .tsv}

\item{suffix}{Suffix to append to all columns [default: NULL]}

\item{suffix.skip}{Number vector indicating which columns to skip when
suffixing [default: 1-2]}
}
\value{
data.frame
}
\description{
Load germline-somatic association statistics from a single cohort
}
