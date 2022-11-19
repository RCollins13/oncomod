#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Constants for RAS modifier project


#' Load RAS Modifier Study Constants
#'
#' Load a subset of constants for RAS modifier analyses
#'
#' @param susbet Vector of constant groups to load See `Details` for options.
#'
#' @details Recognized values for `subset` include:
#' * `colors` : all color palettes used
#' * `scales` : all scales and scale labels
#' * `names` : names of various variables
#' * `all` : load all constants
#'
#' @examples
#' # Load list of color palettes
#' get.constants("colors");
#'
#' # Load scales and names colors
#' get.constants(c("scales", "names"))
#'
#' @export load.constants
#' @export
load.constants <- function(subset){
  # Define colors
  PDAC.colors <- c("dark3" = "#240E29",
                   "dark2" = "#481D53",
                   "dark1" = "#6B2C7D",
                   "main" = "#8F3AA6",
                   "light1" = "#A561B8",
                   "light2" = "#BC89CA",
                   "light3" = "#D2B0DB")
  CRAD.colors <- c("dark3" = "#0E1A29",
                   "dark2" = "#1D3553",
                   "dark1" = "#2C4F7D",
                   "main" = "#3A69A6",
                   "light1" = "#6187B8",
                   "light2" = "#89A5CA",
                   "light3" = "#B0C3DB")
  SKCM.colors <- c("dark3" = "#26290E",
                   "dark2" = "#4B531D",
                   "dark1" = "#717D2C",
                   "main" = "#96A63A",
                   "light1" = "#ABB861",
                   "light2" = "#C0CA89",
                   "light3" = "#D5DBB0")
  male.colors <- c("dark2" = "#2A5869",
                   "dark1" = "#3F839D",
                   "main" = "#54AFD1",
                   "light1" = "#76BFDA",
                   "light2" = "#98CFE3")
  female.colors <- c("dark2" = "#692A4A",
                     "dark1" = "#9D3F6F",
                     "main" = "#D15494",
                     "light1" = "#DA76A9",
                     "light2" = "#E398BF")
  AFR.colors <- c("dark2" = "#796624",
                  "dark1" = "#B59836",
                  "main" = "#F1CB48",
                  "light1" = "#F4D56D",
                  "light2" = "#F7E091")
  AMR.colors <- c("dark2" = "#652223",
                  "dark1" = "#973234",
                  "main" = "#C94345",
                  "light1" = "#D4696A",
                  "light2" = "#DF8E8F")
  EAS.colors <- c("dark2" = "#5B6519",
                  "dark1" = "#889826",
                  "main" = "#B5CA32",
                  "light1" = "#C4D55B",
                  "light2" = "#D3DF84")
  EUR.colors <- c("dark2" = "#4E6774",
                  "dark1" = "#749AAE",
                  "main" = "#9BCDE8",
                  "light1" = "#AFD7ED",
                  "light2" = "#C3E1F1")
  SAS.colors <- c("dark2" = "#4D2D4E",
                  "dark1" = "#734474",
                  "main" = "#995A9B",
                  "light1" = "#AD7BAF",
                  "light2" = "#C29CC3")
  colors <- list(
    "cancer.colors" = c("PDAC" = PDAC.colors[["main"]],
                        "CRAD" = CRAD.colors[["main"]],
                        "SKCM" = SKCM.colors[["main"]]),
    "PDAC.colors" = PDAC.colors,
    "CRAD.colors" = CRAD.colors,
    "SKCM.colors" = SKCM.colors,
    "sex.colors" = c("MALE" = male.colors[["main"]],
                     "FEMALE" = female.colors[["main"]]),
    "MALE.colors" = male.colors,
    "FEMALE.colors" = female.colors,
    "pop.colors" = c("AFR" = AFR.colors[["main"]],
                     "AMR" = AMR.colors[["main"]],
                     "EAS" = EAS.colors[["main"]],
                     "EUR" = EUR.colors[["main"]],
                     "SAS" = SAS.colors[["main"]]),
    "AFR.colors" = AFR.colors,
    "AMR.colors" = AMR.colors,
    "EAS.colors" = EAS.colors,
    "EUR.colors" = EUR.colors,
    "SAS.colors" = SAS.colors,
    "stage.colors" = c("0" = "white",
                       "1" = "#F8FAA7",
                       "2" = "#FFCC66",
                       "3" = "#FE8002",
                       "4" = "#ED3823"))

  # Define scales
  logscale.major <- 10^(-10:10)
  scales <- list(
    "logscale.major" = logscale.major,
    "logscale.major.bp" = 10^(0:9),
    "logscale.major.bp.labels" = c(sapply(c("bp", "kb", "Mb"),
                                          function(suf){paste(c(1, 10, 100), suf, sep="")}),
                                   "1 Gb"),
    "logscale.demi" = as.numeric(sapply(logscale.major, function(e){c(1, 5)*e})),
    "logscale.demi.bp" = as.numeric(sapply(10^(0:9), function(e){c(1, 5)*e})),
    "logscale.demi.bp.labels" = c(paste(c(1, 5, 10, 50, 100, 500), "bp", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "kb", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "Mb", sep=""),
                                  paste(c(1, 5), "Gb", sep="")),
    "logscale.minor" = as.numeric(sapply(logscale.major, function(e){(1:9)*e}))
  )

  # Define names
  all.names <- list(
    "cancer.names.short" = c("PDAC" = "Pancreatic",
                             "CRAD" = "Colorectal",
                             "SKCM" = "Melanoma"),
    "cancer.names.long" = c("PDAC" = "Pancreatic Adenocarcinoma",
                            "CRAD" = "Colorectal Adenocarcinoma",
                            "SKCM" = "Skin Cutaneous Melanoma"),
    "cohort.names.long" = c("TCGA" = "The Cancer Genome Atlas",
                            "PROFILE" = "Dana-Farber Profile"),
    "cohort.names.short" = c("TCGA" = "TCGA",
                             "PROFILE" = "DFCI"),
    "pop.names.short" = c("AFR" = "African",
                          "AMR" = "American",
                          "EAS" = "E. Asian",
                          "EUR" = "European",
                          "SAS" = "S. Asian"),
    "pop.names.long" = c("AFR" = "African/African-American",
                         "AMR" = "Latino/Admixed American",
                         "EAS" = "East Asian",
                         "EUR" = "European",
                         "SAS" = "South Asian"),
    "stage.names" = c("0" = "",
                      "1" = "I",
                      "2" = "II",
                      "3" = "III",
                      "4" = "IV")
  )

  # Assign constants to global environment
  if(length(intersect(subset, c("colors", "all"))) > 0){
    for(variable in names(colors)){
      assign(variable, colors[[variable]], envir=.GlobalEnv)
    }
  }
  if(length(intersect(subset, c("scales", "all"))) > 0){
    for(variable in names(scales)){
      assign(variable, scales[[variable]], envir=.GlobalEnv)
    }
  }
  if(length(intersect(subset, c("names", "all"))) > 0){
    for(variable in names(all.names)){
      assign(variable, all.names[[variable]], envir=.GlobalEnv)
    }
  }
}
