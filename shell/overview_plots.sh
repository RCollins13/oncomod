#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Generate summary plots for cohorts, cancer types, and somatic/germline variants

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export HMFDIR=/data/gusev/USERS/rlc47/HMF
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export CODEDIR=$WRKDIR/../code/oncomod
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in plots plots/overview plots/overview/germline_variants \
              plots/overview/germline_variants/AF_comparisons \
              plots/overview/somatic_variants \
              plots/overview/somatic_variants/freq_comparisons \
              data/plotting; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Ensure most recent version of OncoMod & rCNV2 R packages are installed from source
Rscript -e "install.packages('$CODEDIR/src/OncoModR_0.2.tar.gz', \
                             lib='~/R/x86_64-pc-linux-gnu-library/3.6', \
                             type='source', repos=NULL)"
cd $WRKDIR/../code/rCNV2 && \
git pull && \
cd - && \
Rscript -e "install.packages('$WRKDIR/../code/rCNV2/source/rCNV2_1.0.1.tar.gz', \
                             lib='~/R/x86_64-pc-linux-gnu-library/3.6', \
                             type='source', repos=NULL)"


### Plot patient metadata summaries
$CODEDIR/scripts/plot/plot_pheno_summary.R \
  --cohort-name HMF --metadata $HMFDIR/data/sample_info/HMF.ALL.sample_metadata.tsv.gz \
  --cohort-name DFCI --metadata $PROFILEDIR/data/sample_info/PROFILE.ALL.sample_metadata.tsv.gz \
  --cohort-name TCGA --metadata $TCGADIR/data/sample_info/TCGA.ALL.sample_metadata.tsv.gz \
  --out-prefix $WRKDIR/plots/overview/cohort_summary


### Plot somatic variant summaries
# Collapse all variant frequencies for all cohorts
for cohort in HMF PROFILE TCGA; do
  case $cohort in
    HMF)
      cname=HMF
      ;;
    PROFILE)
      cname=DFCI
      ;;
    TCGA)
      cname=TCGA
      ;;
  esac
  zcat $WRKDIR/data/variant_set_freqs/$cohort.somatic.coding_variants.freq.tsv.gz \
  | sed '1d' | awk -v OFS="\t" -v cohort=$cname '{ print cohort, $0 }'
  zcat $WRKDIR/data/variant_set_freqs/$cohort.somatic.other_variants.freq.tsv.gz \
  | grep -e 'AMP\|DEL' | awk -v OFS="\t" -v cohort=$cname '{ print cohort, $0 }'
done \
| sort -Vk2,2 -k1,1V \
| cat <( zcat $WRKDIR/data/variant_set_freqs/TCGA.somatic.coding_variants.freq.tsv.gz \
         | head -n1 | awk -v OFS="\t" '{ print "cohort", $0 }' ) - \
| gzip -c \
> $TMPDIR/somatic_variant_freqs.tsv.gz
# Build simple table of variant coordinates for each cohort 
for cohort in HMF PROFILE TCGA; do
  case $cohort in
    HMF)
      COHORTDIR=$HMFDIR
      cname=HMF
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      cname=DFCI
      ;;
    TCGA)
      COHORTDIR=$TCGADIR
      cname=TCGA
      ;;
  esac
  bcftools query \
    -f '%ID\t%CHROM\t%POS\n' \
    --regions-file $CODEDIR/refs/RAS_loci.GRCh37.bed.gz \
    $COHORTDIR/data/$cohort.RAS_loci.anno.clean.vcf.gz \
  | awk -v OFS="\t" -v cohort=$cname '{ print cohort, $1, $2, $3 }'
done \
| sort -Vk3,3 -k4,4n -k1,1V | uniq \
| cat <( echo -e "cohort\tvid\tchrom\tpos" ) - \
| gzip -c \
> $TMPDIR/somatic_variant_coords.tsv.gz
# Collapse all variant sets across cohorts
for cohort in HMF PROFILE TCGA; do
  case $cohort in
    HMF)
      cname=HMF
      ;;
    PROFILE)
      cname=DFCI
      ;;
    TCGA)
      cname=TCGA
      ;;
  esac
  zcat $WRKDIR/data/variant_sets/$cohort.somatic.collapsed_coding_csqs.tsv.gz \
  | sed '1d' | awk -v OFS="\t" -v cohort=$cname '{ print cohort, $1, $NF }'
  zcat $WRKDIR/data/variant_sets/$cohort.somatic.other_single_variants.tsv.gz \
  | grep -e 'AMP\|DEL' | awk -v OFS="\t" -v cohort=$cname '{ print cohort, $1, $NF }'
done \
| sort -Vk2,2 -k1,1V | uniq \
| cat <( echo -e "cohort\tset_id\tvids" ) - \
| gzip -c \
> $TMPDIR/variant_set_map.tsv.gz
# Gather necessary plotting data into single file
$CODEDIR/scripts/plot/gather_somatic_ras_data.py \
  --freqs $TMPDIR/somatic_variant_freqs.tsv.gz \
  --variant-coords $TMPDIR/somatic_variant_coords.tsv.gz \
  --variant-set-map $TMPDIR/variant_set_map.tsv.gz \
  --transcript-info $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
  --outfile $WRKDIR/data/plotting/ras_somatic_variants.tsv.gz
# Plot KRAS somatic overview plot
# TODO: implement this
# Scatterplots of inter-cohort somatic frequency correlations for KRAS
$CODEDIR/scripts/plot/plot_somatic_freq_comparisons.R \
  --stats $WRKDIR/data/plotting/ras_somatic_variants.tsv.gz \
  --out-prefix $WRKDIR/plots/overview/somatic_variants/freq_comparisons/somatic_freq_comparisons


### Plot germline variant summaries
# Convert each cohort's AF-annotated VCF to BED without samples
for cohort in TCGA PROFILE HMF; do
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
  svtk vcf2bed -i ALL --no-samples \
    $COHORTDIR/data/$cohort.RAS_loci.anno.clean.wAF.vcf.gz \
    - | bgzip -c \
  > $WRKDIR/data/plotting/$cohort.germline_variants.bed.gz
done
# Plot correlations of germline AFs (inter-cohort & each cohort vs. gnomAD)
$CODEDIR/scripts/plot/plot_germline_AF_comparisons.R \
  --bed $WRKDIR/data/plotting/HMF.germline_variants.bed.gz \
  --name HMF \
  --bed $WRKDIR/data/plotting/PROFILE.germline_variants.bed.gz \
  --name DFCI \
  --bed $WRKDIR/data/plotting/TCGA.germline_variants.bed.gz \
  --name TCGA \
  --out-prefix $WRKDIR/plots/overview/germline_variants/AF_comparisons/AF_comparisons

