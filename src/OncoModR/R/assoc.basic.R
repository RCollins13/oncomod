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
#' @param elig.controls Sample IDs eligible to be reported as 0 \[default: all samples\]
#' @param action Action to apply to query rows; see `Details`
#' @param missing.vid.fill Value to fill for all samples if a specified variant
#' ID from `vids` is not present in `ad` \[default: NA\]
#'
#' @details The `ad` argument accepts either a single data.frame or a list
#' of data.frames. If a list is provided, each matrix will be queried individually
#' and their results will be concatenated.
#'
#' Recognized values for `action` include:
#' * `"verbose"` : return the full query matrix \[default\]
#' * `"any"` : return numeric indicator if any of the query rows are non-zero
#' for each sample
#' * `"all"` : return numeric indicator if all of the query rows are non-zero
#' for each sample
#' * `"sum"` : return the sum of allele dosages for all query rows per sample
#' * `"max"` : return the max of allele dosages for all query rows per sample
#' * `"mean"` : return the mean of allele dosages for all query rows per sample
#'
#' @return numeric vector or data.frame, depending on `action`
#'
#' @seealso [load.ad.matrix], [compress.ad.matrix]
#'
#' @export query.ad.matrix
#' @export
query.ad.matrix <- function(ad, vids, elig.controls=NULL, action="verbose",
                            missing.vid.fill=NA){

  # If ad is a list, call this function recursively on each separately and
  # return the combined query results
  if(inherits(ad, "list")){
    res <- lapply(ad, query.ad.matrix, vids=vids,
                  elig.controls=elig.controls, action=action,
                  missing.vid.fill=missing.vid.fill)
    if(is.data.frame(res[[1]])){
      res <- lapply(res, function(df){df[vids, ]})
      return(as.data.frame(do.call("cbind", res)))
    }else{
      return(unlist(res))
    }
  }

  # Subset to variants of interest, fill missing variants, and return directly if optioned
  sub.df <- ad[which(rownames(ad) %in% vids), ]
  for(vid in setdiff(vids, rownames(sub.df))){
    sub.df[vid, ] <- missing.vid.fill
  }
  if(action == "verbose"){
    return(sub.df)
  }

  # Otherwise, apply various compression strategies prior to returning
  return(compress.ad.matrix(sub.df, action=action, elig.controls=elig.controls))
}


#' Compress allele dosage matrix
#'
#' Compress an allele dosage matrix based on a specified compression function
#'
#' @param ad.df An allele dosage matrix imported by [load.ad.matrix]
#' @param action Action to apply to query rows; see `Details`. Required.
#' @param elig.controls Sample IDs eligible to be reported as 0 \[default: all samples\]
#'
#' @details Recognized values for `action` include:
#' * `"verbose"` : return the full query matrix \[default\]
#' * `"any"` : return numeric indicator if any of the query rows are non-zero
#' for each sample
#' * `"all"` : return numeric indicator if all of the query rows are non-zero
#' for each sample
#' * `"sum"` : return the sum of allele dosages for all query rows per sample
#' * `"max"` : return the max of allele dosages for all query rows per sample
#' * `"mean"` : return the mean of allele dosages for all query rows per sample
#'
#' @return numeric or logical vector, depending on `action`
#'
#' @seealso [load.ad.matrix], [query.ad.matrix]
#'
#' @export compress.ad.matrix
#' @export
compress.ad.matrix <- function(ad.df, action, elig.controls=NULL){
  col.all.na <- apply(ad.df, 2, function(vals){all(is.na(vals))})
  if(action == "any"){
    query.res <- as.numeric(apply(ad.df, 2, function(vals){any(as.logical(vals), na.rm=T)}))
  }else if(action == "all"){
    query.res <- as.numeric(apply(ad.df, 2, function(vals){all(as.logical(vals), na.rm=T)}))
  }else if(action == "sum"){
    query.res <- apply(ad.df, 2, sum, na.rm=T)
  }else if(action == "max"){
    query.res <- apply(ad.df, 2, max, na.rm=T)
  }else if(action == "mean"){
    query.res <- apply(ad.df, 2, mean, na.rm=T)
  }else{
    stop(paste("action ", action, " not recognized as a viable option for compress.ad.matrix", sep="'"))
  }
  query.res[col.all.na] <- NA
  names(query.res) <- colnames(ad.df)
  if(!is.null(elig.controls)){
    query.res[which(query.res == 0 & !(names(query.res) %in% elig.controls))] <- NA
  }
  return(query.res)
}


#' Query genotype quality matrix
#'
#' Extract and compute genotype-conditional mean GQ values per sample
#' from a genotype quality matrix
#'
#' @param gq One or more matrixes of genotype qualities imported by
#' [load.ad.matrix]. See `Details`.
#' @param ad One or more allele dosage matrixes imported by [load.ad.matrix].
#' Order must match `gq`. See `Details`.
#' @param vids Variant IDs to query \[default: use all VIDs in `ad.matrix`\]
#' @param fill.mean Fill missing mean GQs with cohort-wide mean on a
#' genotype-specific basis \[default: TRUE\]
#'
#' @details The `gq` argument accepts either a single data.frame or a list
#' of data.frames. If a list is provided, each matrix will be queried individually
#' and their results will be concatenated.
#'
#' @returns list of two numeric vectors, one for reference genotypes and one for
#' non-reference genotypes
#'
#' @seealso [load.ad.matrix], [query.ad.matrix]
#' @export query.gq.matrix
#' @export
query.gq.matrix <- function(gq, ad, vids=NULL, fill.mean=TRUE){

  # If gq is a list, call this function recursively on each separately and
  # return the combined query results
  if(inherits(gq, "list")){
    if(length(gq) != length(ad)){
      stop("Number of GQ matrixes does not match number of AD matrixes. Stopping.")
    }
    res <- lapply(1:length(gq), function(i){
      query.gq.matrix(gq[[i]], ad[[i]], vids=vids, fill.mean=fill.mean)
    })
    return(list("ref.gq" = unlist(lapply(res, function(l){l$ref.gq})),
                "alt.gq" = unlist(lapply(res, function(l){l$alt.gq}))))
  }

  # Subset AD and GQ matrixes to variants of interest
  if(is.null(vids)){
    vids <- rownames(ad)
  }
  gq <- gq[intersect(rownames(gq), vids), ]
  ad <- ad[intersect(rownames(gq), rownames(ad)), ]

  # If no overlapping variants are found, return NA for all samples
  if(nrow(gq) == 0){
    na.vect <- rep(NA, times=ncol(ad))
    names(na.vect) <- colnames(ad)
    return(list("ref.gq" = na.vect, "alt.gq" = na.vect))
  }

  # Compute mean GQ per sample for reference and non-ref GTs separately
  ref.gq <- sapply(colnames(ad), function(sid){
    mean(gq[which(ad[, sid] == 0), sid], na.rm=T)
  })
  alt.gq <- sapply(colnames(ad), function(sid){
    mean(gq[which(ad[, sid] > 0), sid], na.rm=T)
  })

  # Fill missing mean GQs with cohort-wide means, if optioned
  if(fill.mean){
    ref.gq[which(is.na(ref.gq))] <- mean(ref.gq, na.rm=T)
    alt.gq[which(is.na(alt.gq))] <- mean(alt.gq, na.rm=T)
  }

  return(list("ref.gq" = ref.gq, "alt.gq" = alt.gq))
}


#' Germline-somatic association test
#'
#' Run association test for a single germline-somatic interaction
#'
#' @param y.vals Numeric vector indicating presence or absence of somatic endpoint
#' @param x.vals Numeric vector of germline values to test
#' @param meta Metadata for all samples as loaded by [load.patient.metadata]
#' @param gqs Two-element list of mean ref & alt genotype qualities per sample.
#' If provided, will be included as covariates in association model.
#' @param firth.fallback Attempt to use Firth bias-reduced logistic regression when
#' traditional logistic regression fails to converge or dataset is quasi-separable
#' \[default: TRUE\]
#' @param strict.fallback Implement Firth regression if a standard logit model returns
#' any errors or warnings. Setting this to `FALSE` will only default to Firth
#' regression if logit returns any errors or if the standard error of the genotype
#' coefficient exceeds `nonstrict.se.tolerance` \[default: TRUE\]
#' @param nonstrict.se.tolerance If `strict.fallback` is `FALSE`, only use Firth
#' regression if a standard logit model produces a genotype effect standard error
#' exceeding this value \[default: 10\]
#' @param firth.always Always use Firth regression \[default: FALSE\]
#' @param custom.covariates Character vector of columns from `meta` to include
#' as covariates in the regression model. See `Details`.
#' @param multiPop.min.ac Minimum allele count to include non-European populations
#' @param multiPop.min.freq Minimum variant frequency to include non-European populations
#' @param logistf.maxiter Maximum number of iterations for `logistf()` \[default: 1000\]
#' @param model.suffix String to append to the end of the "model" indicator in
#' the returned vector of results \[default: append nothing\]
#' @param only.return.freqs If set to `TRUE`, no association will be conducted
#' and only frequencies will be reported. All association statistics will be
#' reported as `NA` in the returned vector.
#'
#' @details By default, the following covariates will be included:
#' * Age, sex, and an interaction term for age by sex
#' * Top ten principal components \(PCs\)
#' * Tumor purity
#'
#' Any other covariates to be included must be specified in `custom.covariates`
#'
#' @seealso [load.patient.metadata], [logistf::logistf]
#'
#' @return Vector of test results
#' @export germline.somatic.assoc
#' @export
germline.somatic.assoc <- function(y.vals, x.vals, meta, gqs=NULL,
                                   firth.fallback=TRUE, strict.fallback=TRUE,
                                   nonstrict.se.tolerance=10,
                                   firth.always=FALSE, custom.covariates=c(),
                                   multiPop.min.ac=10,
                                   multiPop.min.freq=0.01,
                                   logistf.maxiter=1000,
                                   model.suffix="",
                                   only.return.freqs=FALSE){
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
  test.df$AGE_BY_SEX <- test.df$AGE_AT_DIAGNOSIS * test.df$SEX

  # Add genotype qualities, if optioned
  if(!is.null(gqs)){
    test.df$REF.GQ <- as.numeric(gqs$ref.gq[rownames(test.df)])
    test.df$ALT.GQ <- as.numeric(gqs$alt.gq[rownames(test.df)])
  }

  # Check if X is allele count or Z-score (for PRS)
  germ.is.ac <- if(all(sapply(x.vals, function(x){as.integer(x) == x}))){TRUE}else{FALSE}

  # Test dataset for quasi- or complete-separability
  # In the case of any zero counts in any of the X by Y matrix,
  # revert to Firth regression
  n.samples <- length(samples)
  somatic.ac <- sum(y.vals, na.rm=T)
  if(!germ.is.ac){
    yes_somatic.germline.ac <- mean(x.vals[which(y.vals > 0)], na.rm=T)
    no_somatic.germline.ac <- mean(x.vals[which(y.vals == 0)], na.rm=T)
    firth <- if(firth.always | somatic.ac < 10){TRUE}else{FALSE}
  }else{
    germline.ac <- sum(x.vals, na.rm=T)
    yes_somatic.germline.ac <- sum(x.vals[which(y.vals > 0)], na.rm=T)
    no_somatic.germline.ac <- sum(x.vals[which(y.vals == 0)], na.rm=T)
    empty.result <- c("samples"=n.samples,
                      "somatic_AC"=somatic.ac,
                      "yes_somatic.germline_AC"=yes_somatic.germline.ac,
                      "no_somatic.germline_AC"=no_somatic.germline.ac,
                      "beta"=NA,
                      "beta_SE"=NA,
                      "z"=NA,
                      "chisq"=NA,
                      "model"=NA,
                      "p"=NA,
                      "EUR_only"=eur.only)
    if(only.return.freqs){
      return(empty.result)
    }
    if(firth.always){
      firth <- TRUE
    }else if(firth.fallback){
      if(length(unique(x.vals)) > 2){
        firth <- if(somatic.ac < 10){TRUE}else{FALSE}
      }else{
        x.by.y <- t(sapply(unique(y.vals), function(y){
          sapply(unique(x.vals), function(x){
            length(which(x.vals[which(y.vals==y)]==x))
          })
        }))
        # Require at least two counts per observed X, Y pair
        # Otherwise, use Firth
        firth <- if(any(x.by.y < 2)){TRUE}else{FALSE}
      }
    }else if(any(c(n.samples, somatic.ac, germline.ac) == 0)){
      return(empty.result)
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

  # Check to ensure number of samples is greater than number of features
  if(nrow(test.df) <= ncol(test.df) - 1){
    return(empty.result)
  }

  # Run association model
  logit.regression <- function(data){
    glm(Y ~ ., data=data, family="binomial")
  }
  firth.regression <- function(data){
    logistf(Y ~ ., data=data, control=logistf.control(maxit=logistf.maxiter), flic=TRUE)
  }
  if(firth){
    fit <- tryCatch(firth.regression(test.df),
                    error=function(e){logit.regression(test.df)})
  }else{
    if(firth.fallback){
      if(strict.fallback){
        fit <- tryCatch(logit.regression(test.df),
                        warning=function(w){firth.regression(test.df)},
                        error=function(e){firth.regression(test.df)})
      }else{
        fit <- tryCatch(logit.regression(test.df),
                        error=function(e){firth.regression(test.df)})
        if(fit$method == "glm.fit"){
          if(summary(fit)$coefficients["X", 2] > nonstrict.se.tolerance){
            fit <- firth.regression(test.df)
          }
        }
      }
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
      model <- paste("flic", model.suffix, sep="")
    }else{
      assoc.res <- as.numeric(summary(fit)$coefficients["X", ])
      z <- assoc.res[3]
      chisq <- NA
      model <- paste("logit", model.suffix, sep="")
    }
  }
  if(assoc.res[2] > nonstrict.se.tolerance){
    stop("Beta SE exceeds permitted SE. Increase nonstrict.se.tolerance to bypass this error")
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

