#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Germline-somatic interaction association analyses

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export CODEDIR=$WRKDIR/../code/ras_modifiers
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Determine LD-independent number of germline variants to test
# Merge VCFs across cohorts
bcftools merge \
  -m none \
  -o $WRKDIR/data/all_cohorts.RAS_loci.vcf.gz \
  -O z \
  $TCGADIR/data/TCGA.RAS_loci.vcf.gz \
  $PROFILEDIR/data/PROFILE.RAS_loci.vcf.gz
# LD prune with PLINK
module load plink/1.90b3
plink \
  --vcf $WRKDIR/data/all_cohorts.RAS_loci.vcf.gz \
  --indep-pairwise 100kb 5 0.2 \
  --recode vcf bgz \
  --out $WRKDIR/data/all_cohorts.RAS_loci.pruned
# Count number of variants retained after LD pruning
wc -l $WRKDIR/data/all_cohorts.RAS_loci.pruned.prune.in
