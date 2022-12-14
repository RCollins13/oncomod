#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Define testable sets of variants for germline-somatic interaction analyses

# Note: intended to be executed on the MGB ERISOne cluster


##################
### Basic prep ###
##################

### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export CODEDIR=$WRKDIR/../code/ras_modifiers
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/variant_sets data/variant_set_freqs; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Generate maps of somatic variants with the same protein consequence
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
  esac
  for context in germline somatic; do
    case $context in
      germline)
        subset="RAS_loci"
        ;;
      somatic)
        subset="somatic_variants"
        ;;
    esac
    bsub -q normal \
      -J collapse_coding_csqs_${cohort}_$context \
      -o $WRKDIR/LSF/logs/collapse_coding_csqs_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/collapse_coding_csqs_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/vep2csqTable.py \
        $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz \
        $WRKDIR/data/variant_sets/$cohort.$context.collapsed_coding_csqs.tsv.gz"
  done
done


### Generate somatic and germline variant sets for burden testing
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
  esac
  for context in somatic germline; do
    # bsub -q normal \
    #   -J generate_variant_sets_${cohort}_$subset \
    #   -o $WRKDIR/LSF/logs/generate_variant_sets_${cohort}_$subset.log \
    #   -e $WRKDIR/LSF/logs/generate_variant_sets_${cohort}_$subset.err \
    #   TBD
  done
done
#DEV:
# spec for .json input to make_variant_sets.py
# {
#   // Each set gets a unique identifier
#   set_id:
#   // Within each set is an array of sets of criteria.
#   // The relationship between sets of criteria is assumed to be AND within 
#   // each object (curly brackets) and assumed to be OR between objects
#   [
#     // Each key : value pair specifies a single criterion to require
#     // The value is provided as an array of value and equality
#     // Value can be numeric, boolean, or string
#     // Equality should be two-letter string shorthand for comparison to apply
#     // examples: eq, gt, lt, ge, le
#     {
#       key : [value, equality],
#       key : [value, equality]
#     },
#     {
#       key : [value, equality],
#       key : [value, equality]
#     }
#   ]
# }


### Compute AC and AF matrixes by cancer type for all variants & variant sets in each cohort
# 1. Coding variants by protein consequence
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
  esac
  for context in germline somatic; do
    case $context in
      germline)
        subset="RAS_loci"
        max_an=2
        ;;
      somatic)
        subset="somatic_variants"
        max_an=1
        ;;
    esac
    bsub -q normal -J get_coding_freqs_${cohort}_$context \
      -o $WRKDIR/LSF/logs/get_coding_freqs_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/get_coding_freqs_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/calc_variant_set_freqs.py \
         --sets-tsv $WRKDIR/data/variant_sets/$cohort.$context.collapsed_coding_csqs.tsv.gz \
         --dosage-tsv $COHORTDIR/data/$cohort.$subset.dosage.tsv.gz \
         --sample-metadata $COHORTDIR/data/sample_info/TCGA.ALL.sample_metadata.tsv.gz \
         --max-an $max_an \
         --outfile $WRKDIR/data/variant_set_freqs/$cohort.$context.coding_variants.freq.tsv.gz"
  done
done
# 2. All other single variants after excluding coding variants above
# 3. All variant sets
# TODO: implement this
