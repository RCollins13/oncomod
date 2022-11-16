#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Constants for RAS modifier project


#' Load RAS modifier study constants
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
  PDAC.colors <- list("dark3" = "#220829",
                      "dark2" = "#450F53",
                      "dark1" = "#68167D",
                      "main" = "#8A1EA6",
                      "light1" = "#A14BB8",
                      "light2" = "#B978CA",
                      "light3" = "#D0A5DB")
  CRAD.colors <- list("dark3" = "#081629",
                      "dark2" = "#0F2D53",
                      "dark1" = "#16437D",
                      "main" = "#1E59A6",
                      "light1" = "#4B7AB8",
                      "light2" = "#789BCA",
                      "light3" = "#A5BDDB")
  SKCM.colors <- list("dark3" = "#252908",
                      "dark2" = "#5D6719",
                      "dark1" = "#6E7D16",
                      "main" = "#92A61E",
                      "light1" = "#A8B84B",
                      "light2" = "#BECA78",
                      "light3" = "#D3DBA5")
  male.colors <- list("dark2" = "#0F2D53",
                      "dark1" = "#2C7E9D",
                      "main" = "#3BA8D1",
                      "light1" = "#62B9DA",
                      "light2" = "#89CBE3")
  female.colors <- list("dark2" = "#691D44",
                        "dark1" = "#9D2C66",
                        "main" = "#D13A88",
                        "light1" = "#DA61A0",
                        "light2" = "#E389B8")
  AFR.colors <- list("dark2" = "#796624",
                     "dark1" = "#B59836",
                     "main" = "#F1CB48",
                     "light1" = "#F4D56D",
                     "light2" = "#F7E091")
  AMR.colors <- list("dark2" = "#652223",
                     "dark1" = "#973234",
                     "main" = "#C94345",
                     "light1" = "#D4696A",
                     "light2" = "#DF8E8F")
  EAS.colors <- list("dark2" = "#5B6519",
                     "dark1" = "#889826",
                     "main" = "#B5CA32",
                     "light1" = "#C4D55B",
                     "light2" = "#D3DF84")
  EUR.colors <- list("dark2" = "#4E6774",
                     "dark1" = "#749AAE",
                     "main" = "#9BCDE8",
                     "light1" = "#AFD7ED",
                     "light2" = "#C3E1F1")
  SAS.colors <- list("dark2" = "#4D2D4E",
                     "dark1" = "#734474",
                     "main" = "#995A9B",
                     "light1" = "#AD7BAF",
                     "light2" = "#C29CC3")
  colors <- list("cancer.colors" = list("PDAC" = PDAC.colors$main,
                                        "CRAD" = CRAD.colors$main,
                                        "SKCM" = SKCM.colors$main),
                 "PDAC.colors" = PDAC.colors,
                 "CRAD.colors" = CRAD.colors,
                 "SKCM.colors" = SKCM.colors,
                 "sex.colors" = list("male" = male.colors$main,
                                     "female" = female.colors$main),
                 "male.colors" = male.colors,
                 "female.colors" = female.colors,
                 "pop.colors" = list("AFR" = AFR.colors$main,
                                     "AMR" = AMR.colors$main,
                                     "EAS" = EAS.colors$main,
                                     "EUR" = EUR.colors$main,
                                     "SAS" = SAS.colors$main),
                 "AFR.colors" = AFR.colors,
                 "AMR.colors" = AMR.colors,
                 "EAS.colors" = EAS.colors,
                 "EUR.colors" = EUR.colors,
                 "SAS.colors" = SAS.colors)

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
  rasmod.names <- list(
    "cancer.names.short" = c("PDAC" = "Pancreatic",
                             "CRAD" = "Colorectal",
                             "SKCM" = "Melanoma"),
    "cancer.names.long" = c("PDAC" = "Pancreatic Adenocarcinoma",
                            "CRAD" = "Colorectal Adenocarcinoma",
                            "SKCM" = "Skin Cutaneous Melanoma"),
    "cohort.names" = c("TCGA" = "The Cancer Genome Atlas",
                       "PROFILE" = "Dana-Farber Profile")
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
    for(variable in names(rasmod.names)){
      assign(variable, rasmod.names[[variable]], envir=.GlobalEnv)
    }
  }
}
