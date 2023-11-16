#!/usr/bin/env Rscript

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Script to extract all Tier 1 KRAS mutations reported for a subset of PROFILE samples


# Helper function to sort list of protein mutations by residue number and allele
sortMuts <- function(mut.list){
  residue <- as.numeric(gsub("[A-z.]", "", mut.list))
  aa <- sapply(mut.list, function(m){unlist(strsplit(m, split="[0-9]+"))[2]})
  mdf <- data.frame("r"=residue, "a"=aa)
  mut.list[with(mdf, order(r, a))]
}


###########
# RScript #
###########
# Parse command line arguments and options
require(argparse, quietly=TRUE)
options(stringsAsFactors=F, scipen=1000)
parser <- ArgumentParser(description=paste("Gather all tier 1 pathogenic KRAS",
                                           "mutations reported for a subset of",
                                           "PROFILE samples"))
parser$add_argument("--mutations-csv", metavar=".csv", type="character",
                    help="OncDRS genomic mutation results", required=TRUE)
parser$add_argument("--sample-ids", metavar=".txt", type="character",
                    help="list of samples to be included")
parser$add_argument("--id-map", metavar=".tsv", type="character",
                    help=paste(".tsv mapping various forms of PROFILE IDs.",
                               "Required if --sample-ids is also provided."))
parser$add_argument("--outfile", metavar=".txt", type="character", required=TRUE,
                    help="output .txt file for list of all tier 1 mutations")
args <- parser$parse_args()

# Load mutation data
muts <- read.table(args$mutations_csv, sep=",", header=T, comment.char="")

# Subset to tier 1 KRAS mutations
muts <- muts[which(muts$PATHOLOGIST_TIER == 1
                   & muts$CANONICAL_GENE == "KRAS"), ]

# Subset to patients of interest
if(!is.null(args$sample_ids)){
  if(is.null(args$id_map)){
    stop("Error: must supply --id-map if --sample-ids is specified")
  }else{
    idmap <- read.table(args$id_map, sep="\t", comment.char="", header=T)
  }
  sids <- read.table(args$sample_ids, header=F)[, 1]
  blids <- idmap$SAMPLE_ACCESSION_NBR[which(idmap$PBP %in% sids)]
  muts <- muts[which(muts$SAMPLE_ACCESSION_NBR %in% blids), ]
}

# Get list of all unique protein consequences, sort by residue number, and write to file
tier1 <- unique(muts$CANONICAL_PROTEIN_CHANGE)
write.table(sortMuts(tier1), args$outfile, row.names=F, col.names=F, sep="\t", quote=F)
