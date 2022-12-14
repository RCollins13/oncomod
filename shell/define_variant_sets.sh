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
tabix -H $PROFILEDIR/data/PROFILE.RAS_loci.anno.clean.vcf.gz > $TMPDIR/test.vcf
zcat $PROFILEDIR/data/PROFILE.RAS_loci.anno.clean.vcf.gz | grep -ve '^#' | grep -e 'enhancer\|promoter\|SpliceAI\|GTEx' >> $TMPDIR/test.vcf
bgzip -f $TMPDIR/test.vcf
tabix -f $TMPDIR/test.vcf.gz
$TMPDIR/populate_variant_sets.py \
  --vcf $TMPDIR/test.vcf.gz \
  --sets-json $TMPDIR/variant_set_criteria.germline.json \
  --outfile $TMPDIR/test.burden_sets.tsv


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
    for suf in err log; do
      logfile=$WRKDIR/LSF/logs/get_coding_freqs_${cohort}_$context.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
      -J get_coding_freqs_${cohort}_$context \
      -o $WRKDIR/LSF/logs/get_coding_freqs_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/get_coding_freqs_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/calc_variant_set_freqs.py \
         --sets-tsv $WRKDIR/data/variant_sets/$cohort.$context.collapsed_coding_csqs.tsv.gz \
         --dosage-tsv $COHORTDIR/data/$cohort.$subset.dosage.tsv.gz \
         --sample-metadata $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
         --max-an $max_an \
         --outfile $WRKDIR/data/variant_set_freqs/$cohort.$context.coding_variants.freq.tsv.gz"
  done
done
# 2. All other single variants after excluding coding variants above
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
    # First, need to get list of variant IDs not already included in coding sets above
    zcat $WRKDIR/data/variant_sets/$cohort.$context.collapsed_coding_csqs.tsv.gz \
    | sed '1d' | awk -v FS="\t" '{ print $NF }' | sed 's/,/\n/g' | sort | uniq \
    > $TMPDIR/$cohort.$context.coding_vids.list
    bcftools query \
      --exclude "ID = @$TMPDIR/$cohort.$context.coding_vids.list" \
      --format '%ID\n' \
      $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz \
    | sort -V | uniq | awk -v OFS="\t" '{ print $1, $1 }' \
    | cat <( echo -e "set_id\tvids" ) - | gzip -c \
    > $WRKDIR/data/variant_sets/$cohort.$context.other_single_variants.tsv.gz
    # After this list is compiled, now we can submit the frequency collection task
    for suf in err log; do
      logfile=$WRKDIR/LSF/logs/get_other_freqs_${cohort}_$context.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
      -J get_other_freqs_${cohort}_$context \
      -o $WRKDIR/LSF/logs/get_other_freqs_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/get_other_freqs_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/calc_variant_set_freqs.py \
         --sets-tsv $WRKDIR/data/variant_sets/$cohort.$context.other_single_variants.tsv.gz \
         --dosage-tsv $COHORTDIR/data/$cohort.$subset.dosage.tsv.gz \
         --sample-metadata $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
         --max-an $max_an \
         --outfile $WRKDIR/data/variant_set_freqs/$cohort.$context.other_variants.freq.tsv.gz"
  done
done
# 3. All variant sets
# TODO: implement this
