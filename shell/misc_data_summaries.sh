#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Helper code to gather useful tables and data snippets for slides


##################
### Basic prep ###
##################

### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export HMFDIR=/data/gusev/USERS/rlc47/HMF
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export CODEDIR=$WRKDIR/../code/oncomod
cd $WRKDIR


#########################################################
### Get count of germline SNVs/indels in various loci ###
#########################################################

# KRAS locus by cohort
for cohort in PROFILE HMF TCGA; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
    HMF)
      COHORTDIR=$HMFDIR
      ;;
  esac

  echo $cohort
  
  for cancer in PDAC CRAD LUAD; do
    bcftools query \
      -f '%CHROM\_%POS\_%REF\_%ALT\n' \
      --include "INFO/${cancer}_AC > 0" \
      --regions-file <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep KRAS ) \
      $COHORTDIR/data/$cohort.RAS_loci.anno.clean.wAF.vcf.gz \
    > $TMPDIR/$cohort.$cancer.KRAS_cis_alleles.vids.list
    cat $TMPDIR/$cohort.$cancer.KRAS_cis_alleles.vids.list | wc -l
  done

  bcftools query \
    -f '%CHROM\_%POS\_%REF\_%ALT\n' \
      --include "AC > 0" \
      --regions-file <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep KRAS ) \
      $COHORTDIR/data/$cohort.RAS_loci.anno.clean.wAF.vcf.gz \
  > $TMPDIR/$cohort.ALL.KRAS_cis_alleles.vids.list
  cat $TMPDIR/$cohort.ALL.KRAS_cis_alleles.vids.list | wc -l
done | paste - - - - -

# KRAS locus, intersection of all cohorts
# Note: must have run the code block directly above
for cancer in PDAC CRAD LUAD ALL; do
  for cohort in PROFILE HMF TCGA; do
    cat $TMPDIR/$cohort.$cancer.KRAS_cis_alleles.vids.list
  done | sort -V | uniq -c | awk '{ if ($1==3) print }' | wc -l
done | paste <( echo "Intersection" ) - - - -

# KRAS locus, union of all cohorts
# Note: must have run the code block directly above
for cancer in PDAC CRAD LUAD ALL; do
  for cohort in PROFILE HMF TCGA; do
    cat $TMPDIR/$cohort.$cancer.KRAS_cis_alleles.vids.list
  done | sort -V | uniq | wc -l
done | paste <( echo "Union" ) - - - -


##################################################################
### Get count of eligible "control" patients per cohort/cancer ###
##################################################################
for cancer in PDAC CRAD LUAD; do
  echo $cancer
  for cohort in PROFILE HMF TCGA; do
    case $cohort in
      TCGA)
        COHORTDIR=$TCGADIR
        sample_field="donors"
        ;;
      PROFILE)
        COHORTDIR=$PROFILEDIR
        sample_field="samples"
        ;;
      HMF)
        COHORTDIR=$HMFDIR
        sample_field="samples"
        ;;
    esac
    elig_samps=$COHORTDIR/data/sample_info/$cohort.$cancer.$sample_field.list
    fgrep -wf $elig_samps \
      $COHORTDIR/data/sample_info/$cohort.ALL.eligible_controls.list \
    | wc -l | addcom
  done | paste -s -
done | paste - -

