#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Ensure non-standard dependencies are available when package is loaded


.onLoad <- function(libname, pkgname){

  # Set useful global constants
  options(scipen=1000, stringsAsFactors=F, family="sans")

  # Check to make sure RLCtools is available
  if(!require(RLCtools)){
    stop(paste("Dependency `RLCtools` is required for `OncoModR` but is not found.\n",
               "For more info, see https://github.com/RCollins13/RLCtools\n", sep=""))
  }
}
