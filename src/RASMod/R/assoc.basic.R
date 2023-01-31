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
#' @param ad One or more allele dosage matrix(es) as imported by
#' [load.ad.matrix]. See `Details`
#' @param vids Variant IDs to query
#' @param action Action to apply to query rows; see `Details`
#'
#' @details The `ad` argument accepts either a single data.frame or a list
#' of data.frames. If a list is provided, each matrix will be queried individually
#' and their results will be concatenated.
#'
#' Recognized values for `action` include:
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

  # If ad is a list, call this function recursively on each separately and
  # return the combined query results
  if(inherits(ad, "list")){
    res <- lapply(ad, query.ad.matrix, vids=vids, action=action)
    if(is.data.frame(res[[1]])){
      return(as.data.frame(do.call("rbind", res)))
    }else{
      return(unlist(res))
    }
  }

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
    query.res <- apply(sub.df, 2, sum, na.rm=T)
  }
  query.res[col.all.na] <- NA
  names(query.res) <- colnames(sub.df)
  return(query.res)
}


#' Germline-somatic association test
#'
#' Run association test for a single germline-somatic interaction
#'
#' @param y.vals Numeric vector indicating presence or absence of somatic endpoint
#' @param x.vals Numeric vector of germline values to test
#' @param samples Vector of sample IDs to consider
#' @param meta Metadata for all samples as loaded by [load.patient.metadata]
#' @param firth.fallback Attempt to use Firth bias-reduced logistic regression when
#' traditional logistic regression fails to converge or dataset is quasi-separable
#' \[default: TRUE\]
#' @param firth.always Always use Firth regression \[default: FALSE\]
#' @param custom.covariates Character vector of columns from `meta` to include
#' as covariates in the regression model. See `Details`.
#' @param multiPop.min.ac Minimum allele count to include non-European populations
#' @param multiPop.min.freq Minimum variant frequency to include non-European populations
#'
#' @details By default, the following covariates will be included:
#' * Age, sex, and an interaction term for age by sex
#' * Top ten principal components \(PCs\)
#' * Tumor purity
#'
#' Any other covariates to be included must be specified in `custom.covariates`
#'
#' @seealso [load.patient.metadata]
#'
#' @return Vector of test results
#' @export germline.somatic.assoc
#' @export
germline.somatic.assoc <- function(y.vals, x.vals, samples, meta,
                                   firth.fallback=TRUE, firth.always=FALSE,
                                   custom.covariates=c(), multiPop.min.ac=10,
                                   multiPop.min.freq=0.01){
  # Ensure Firth package is loaded
  require(logistf, quietly=TRUE)

  # Enforce strict intersection of non-NA samples between germline & somatic
  som.samples <- names(y.vals)[which(!is.na(y.vals))]
  germ.samples <- names(x.vals)[which(!is.na(x.vals))]
  samples <- intersect(som.samples, germ.samples)
  y.vals <- y.vals[samples]
  x.vals <- x.vals[samples]

  # If germline or somatic variants are <1% freq / <10 counts, restrict to
  # European samples to protect against pop strat
  eur.only <- FALSE
  if(sum(x.vals) < multiPop.min.ac
     | sum(x.vals) / length(samples) < multiPop.min.freq
     | sum(y.vals) < multiPop.min.ac
     | sum(y.vals) / length(samples) < multiPop.min.freq){
    eur.samples <- rownames(meta[which(meta$POPULATION == "EUR"), ])
    samples <- intersect(samples, eur.samples)
    x.vals <- x.vals[samples]
    y.vals <- y.vals[samples]
    eur.only <- TRUE
  }

  # Ensure non-zero counts for X and Y
  if(length(samples) == 0
     | length(table(x.vals)) < 2
     | length(table(y.vals)) < 2){
    return(NULL)
  }

  # Construct test df from y, x, and meta
  cov.to.keep <- Reduce(union,
                        list(c("AGE_AT_DIAGNOSIS", "SEX", "TUMOR_PURITY"),
                            custom.covariates,
                            colnames(meta)[grep("^PC[0-9]", colnames(meta))]))
  test.df <- meta[samples, intersect(cov.to.keep, colnames(meta))]
  test.df$SEX <- as.numeric(test.df$SEX == "MALE")
  test.df <- cbind(data.frame("Y" = y.vals, "X" = x.vals, row.names=names(y.vals)),
                   test.df)

  # Check if X is allele count or Z-score (for PRS)
  germ.is.ac <- if(all(is.integer(x.vals))){TRUE}else{FALSE}

  # Test dataset for quasi- or complete-separability
  # In the case of any zero counts in any of the X by Y matrix,
  # revert to Firth regression
  n.samples <- length(samples)
  somatic.ac <- sum(y.vals, na.rm=T)
  if(!germ.is.ac){
    yes_somatic.germline.ac <- mean(x.vals[which(y.vals > 0)], na.rm=T)
    no_somatic.germline.ac <- mean(x.vals[which(y.vals == 0)], na.rm=T)
    firth <- if(somatic.ac < 10){TRUE}else{FALSE}
  }else{
    germline.ac <- sum(x.vals, na.rm=T)
    yes_somatic.germline.ac <- sum(x.vals[which(y.vals > 0)], na.rm=T)
    no_somatic.germline.ac <- sum(x.vals[which(y.vals == 0)], na.rm=T)
    if(firth.always){
      firth <- TRUE
    }else if(firth.fallback){
      x.by.y <- t(sapply(unique(y.vals), function(y){
        sapply(unique(x.vals), function(x){
          length(which(x.vals[which(y.vals==y)]==x))
        })
      }))
      # Require at least two counts per observed X, Y pair
      # Otherwise, use Firth
      if(any(x.by.y < 2)){
        firth <- TRUE
      }else{
        firth <- FALSE
      }
    }else if(any(c(n.samples, somatic.ac, germline.ac) == 0)){
      return(c("samples"=n.samples,
               "somatic_AC"=somatic.ac,
               "yes_somatic.germline_AC"=yes_somatic.germline.ac,
               "no_somatic.germline_AC"=no_somatic.germline.ac,
               "beta"=NA,
               "beta_SE"=NA,
               "z"=NA,
               "chisq"=NA,
               "model"=NA,
               "p"=NA,
               "EUR_only"=eur.only))
    }else{
      firth <- FALSE
    }
  }

  # Drop covariates with no informative observations
  useless.cov <- which(apply(test.df, 2, function(vals){
    all(is.na(vals)) | (length(table(vals)) < 2)
    }))
  if(length(useless.cov > 0)){
    test.df <- test.df[, -useless.cov]
  }

  # Standard normalize all covariates
  test.df[, -c(1:2)] <- apply(test.df[, -c(1:2)], 2, scale)

  # Run association model
  logit.regression <- function(data){
    glm(Y ~ . + (AGE_AT_DIAGNOSIS * SEX), data=data, family="binomial")
  }
  firth.regression <- function(data){
    logistf(Y ~ . + (AGE_AT_DIAGNOSIS * SEX), data=data,
            control=logistf.control(maxit=100), flic=TRUE)
  }
  if(firth){
    fit <- tryCatch(firth.regression(test.df),
                    error=function(e){logit.regression(test.df)})
  }else{
    if(firth.fallback){
      fit <- tryCatch(logit.regression(test.df),
                      warning=function(w){firth.regression(test.df)},
                      error=function(e){firth.regression(test.df)})
    }else{
      fit <- logit.regression(test.df)
    }
  }
  if(fit$method != "glm.fit"){
    firth <- TRUE
  }

  # Extract association stats for germline variants
  if(is.na(fit$coefficients["X"])){
    assoc.res <- rep(NA, 4)
    model <- chisq <- z <- NA
  }else{
    if(firth){
      assoc.res <- as.numeric(c(fit$coefficients["X"],
                                sqrt(diag(vcov(fit)))["X"],
                                qchisq(1-fit$prob, df=1)["X"],
                                fit$prob["X"]))
      chisq <- assoc.res[3]
      z <- NA
      model <- "flic"
    }else{
      assoc.res <- as.numeric(summary(fit)$coefficients["X", ])
      z <- assoc.res[3]
      chisq <- NA
      model <- "logit"
    }
  }
  # Return summary vector
  res <- c("samples"=n.samples,
           "somatic_AC"=somatic.ac,
           "yes_somatic.germline_AC"=yes_somatic.germline.ac,
           "no_somatic.germline_AC"=no_somatic.germline.ac,
           "beta"=assoc.res[1],
           "beta_SE"=assoc.res[2],
           "z"=z,
           "chisq"=chisq,
           "model"=model,
           "p"=assoc.res[4],
           "EUR_only"=eur.only)
  return(res)
}
