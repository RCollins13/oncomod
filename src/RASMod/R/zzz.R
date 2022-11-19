#!/usr/bin/env R

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Actions to perform upon loading the RASMod package

.onLoad <- function(libname, pkgname){
  options(scipen=1000, stringsAsFactors=F, family="sans")
}

# .onAttach <- function(libname, pkgname){
#   load.constants("all")
# }
