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
for SUBDIR in data data/variant_set_freqs/filtered; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Determine LD-independent number of germline variants to test
# Subset VCFs per cohort to non-rare (AC≥10) variants
bcftools view \
  --min-ac 10 \
  -o $TCGADIR/data/TCGA.RAS_loci.ac10plus.vcf.gz \
  -O z \
  $TCGADIR/data/TCGA.RAS_loci.vcf.gz
tabix -p vcf -f $TCGADIR/data/TCGA.RAS_loci.ac10plus.vcf.gz
bcftools view \
  --min-ac 10 \
  -o $PROFILEDIR/data/PROFILE.RAS_loci.ac10plus.vcf.gz \
  -O z \
  $PROFILEDIR/data/PROFILE.RAS_loci.vcf.gz
tabix -p vcf -f $PROFILEDIR/data/PROFILE.RAS_loci.ac10plus.vcf.gz
# Merge AC≥10 VCFs across cohorts
bcftools merge \
  -m none \
  -o $WRKDIR/data/all_cohorts.RAS_loci.ac10plus.vcf.gz \
  -O z \
  $TCGADIR/data/TCGA.RAS_loci.ac10plus.vcf.gz \
  $PROFILEDIR/data/PROFILE.RAS_loci.ac10plus.vcf.gz
# LD prune with PLINK
module load plink/1.90b3
plink \
  --vcf $WRKDIR/data/all_cohorts.RAS_loci.ac10plus.vcf.gz \
  --indep-pairwise 100kb 5 0.2 \
  --recode vcf bgz \
  --out $WRKDIR/data/all_cohorts.RAS_loci.pruned
# Count number of variants retained after LD pruning
wc -l $WRKDIR/data/all_cohorts.RAS_loci.pruned.prune.in


### Filter somatic variants to define list of conditions to test
# 1. Frequent RAS mutations
for cohort in TCGA PROFILE; do
  # Coding (collapsed by consequence)
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.somatic.coding_variants.freq.tsv.gz
  zcat $WRKDIR/data/variant_sets/$cohort.somatic.collapsed_coding_csqs.tsv.gz \
  | grep -e 'KRAS\|NRAS\|HRAS' | cut -f1 \
  | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-freq 0.01 \
    --outfile $WRKDIR/data/variant_set_freqs/$cohort.somatic.${context}_variants.freq.1pct.tsv.gz
  # Noncoding (individual variants)
      

  done
done
# 2. Recurrently mutated RAS codons
# TODO: implement this
# 3. Functional mutation sets
# TODO: implement this
# 4. Frequent RAS co-mutation pairs
# TODO: implement this
# 5. RAS + other gene co-mutation pairs
# TODO: implement this
# 6. RAS + RAS signaling co-mutation pairs
