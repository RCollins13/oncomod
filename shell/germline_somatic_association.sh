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
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.coding_variants.freq.1pct.tsv.gz

  # Noncoding (individual variants)
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
  esac
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.somatic.other_variants.freq.tsv.gz
  tabix \
    -R $WRKDIR/../refs/RAS_genes.bed.gz \
    $COHORTDIR/data/$cohort.somatic_variants.anno.clean.vcf.gz \
  | cut -f3 | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-freq 0.01 \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.other_variants.freq.1pct.tsv.gz
done

# 2. Recurrently mutated RAS codons
for cohort in TCGA PROFILE; do
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.somatic.recurrently_mutated_codons.freq.tsv.gz
  zcat $WRKDIR/data/variant_sets/$cohort.somatic.recurrently_mutated_codons.tsv.gz \
  | fgrep -wf \
    <( zcat $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
       | awk -v FS="\t" '{ if ($3 == "NRAS" || $3 == "HRAS" || $3 == "KRAS") print $1 }' \
       | cut -f1 -d\. ) \
  | cut -f1 | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-freq 0.01 \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.recurrently_mutated_codons.freq.1pct.tsv.gz
done

# 3. Functional mutation sets
for cohort in TCGA PROFILE; do
  zcat $WRKDIR/data/variant_set_freqs/$cohort.somatic.burden_sets.freq.tsv.gz \
  | grep -e '^set_id\|^NRAS_\|^HRAS_\|^KRAS_' \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-freq 0.01 \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.burden_sets.freq.1pct.tsv.gz
done

# 4. Frequent co-mutation pairs involving RAS
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
  esac
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.somatic.gene_comutations.freq.tsv.gz
  tabix \
    -R $WRKDIR/../refs/RAS_genes.bed.gz \
    $COHORTDIR/data/$cohort.somatic_variants.anno.clean.vcf.gz \
  | cut -f3 | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-freq 0.01 \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.comutations.freq.1pct.tsv.gz
done

# 5. Frequent RAS intra-gene (RAS+other) co-mutation pairs involving burden sets in other genes
# TODO: implement this
# 6. RAS + RAS signaling co-mutation pairs


### Summarize somatic conditions to test as endpoints for association
# $CODEDIR/scripts/germline_somatic_assoc/summarize_somatic_endpoints.py \
$TMPDIR/summarize_somatic_endpoints.py \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.coding_variants.freq.1pct.tsv.gz \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.coding_variants.freq.1pct.tsv.gz \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.other_variants.freq.1pct.tsv.gz \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.other_variants.freq.1pct.tsv.gz \
  --codons $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.recurrently_mutated_codons.freq.1pct.tsv.gz \
  --codons $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.recurrently_mutated_codons.freq.1pct.tsv.gz \
  --comutations $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.comutations.freq.1pct.tsv.gz \
  --comutations $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.comutations.freq.1pct.tsv.gz \
  --transcript-info $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz
