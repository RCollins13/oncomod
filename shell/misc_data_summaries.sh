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

# Variants in RAS-driven cancer-associated GWAS loci by cohort
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
      --regions-file <( zcat $CODEDIR/refs/GWAS_catalog_loci.CRAD_PDAC_LUAD.GRCh37.bed.gz ) \
      $COHORTDIR/data/$cohort.RAS_loci.anno.clean.wAF.vcf.gz \
    > $TMPDIR/$cohort.$cancer.GWAS_loci.vids.list
    cat $TMPDIR/$cohort.$cancer.GWAS_loci.vids.list | wc -l
  done

  bcftools query \
    -f '%CHROM\_%POS\_%REF\_%ALT\n' \
      --include "AC > 0" \
      --regions-file <( zcat $CODEDIR/refs/GWAS_catalog_loci.CRAD_PDAC_LUAD.GRCh37.bed.gz ) \
      $COHORTDIR/data/$cohort.RAS_loci.anno.clean.wAF.vcf.gz \
  > $TMPDIR/$cohort.ALL.GWAS_loci.vids.list
  cat $TMPDIR/$cohort.ALL.GWAS_loci.vids.list | wc -l
done | paste - - - - -

# RAS-driven cancer-associated GWAS loci, intersection of all cohorts
# Note: must have run the code block directly above
for cancer in PDAC CRAD LUAD ALL; do
  for cohort in PROFILE HMF TCGA; do
    cat $TMPDIR/$cohort.$cancer.GWAS_loci.vids.list
  done | sort -V | uniq -c | awk '{ if ($1==3) print }' | wc -l
done | paste <( echo "Intersection" ) - - - -

# RAS-driven cancer-associated GWAS loci, union of all cohorts
# Note: must have run the code block directly above
for cancer in PDAC CRAD LUAD ALL; do
  for cohort in PROFILE HMF TCGA; do
    cat $TMPDIR/$cohort.$cancer.GWAS_loci.vids.list
  done | sort -V | uniq | wc -l
done | paste <( echo "Union" ) - - - -

# All KRAS pathway genes, any coding consequence
# 1. Collect variant IDs corresponding to nonsynonymous variants in each gene
#    KRAS locus by cohort and cancer type
for cohort in PROFILE HMF TCGA; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      sample_field=donors
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      sample_field=samples
      ;;
    HMF)
      COHORTDIR=$HMFDIR
      sample_field=samples
      ;;
  esac
  if [ -e $COHORTDIR/data/ras_pathway_coding_vid_lists ]; then
    rm -rf $COHORTDIR/data/ras_pathway_coding_vid_lists
  fi
  mkdir $COHORTDIR/data/ras_pathway_coding_vid_lists
  for cancer in PDAC CRAD LUAD ALL; do
    while read gene; do
      cat << EOF > $COHORTDIR/LSF/scripts/$cohort.$cancer.$gene.get_nonsyn_vids.sh
bcftools view \
  --samples-file $COHORTDIR/data/sample_info/$cohort.$cancer.$sample_field.list \
  --min-ac 1 \
  $COHORTDIR/data/$cohort.RAS_loci.anno.clean.vcf.gz \
| bcftools +split-vep - \
  -f '%CHROM\_%POS\_%REF\_%ALT\n' \
  -c SYMBOL \
  -s worst:incomplete_terminal_codon+ \
  -i "SYMBOL == \"$gene\"" \
> $COHORTDIR/data/ras_pathway_coding_vid_lists/$cohort.$cancer.$gene.nonsyn.vids.list
EOF
      chmod a+x $COHORTDIR/LSF/scripts/$cohort.$cancer.$gene.get_nonsyn_vids.sh
      for suf in err log; do
        if [ -e $COHORTDIR/LSF/logs/$cohort.$cancer.$gene.get_nonsyn_vids.$suf ]; then
          rm $COHORTDIR/LSF/logs/$cohort.$cancer.$gene.get_nonsyn_vids.$suf
        fi
      done
      bsub -q normal -sla miket_sc \
        -J cohort.$cancer.$gene.get_nonsyn_vids \
        -o $COHORTDIR/LSF/logs/$cohort.$cancer.$gene.get_nonsyn_vids.log \
        -e $COHORTDIR/LSF/logs/$cohort.$cancer.$gene.get_nonsyn_vids.err \
        $COHORTDIR/LSF/scripts/$cohort.$cancer.$gene.get_nonsyn_vids.sh
    done < $CODEDIR/refs/NCI_RAS_pathway.genes.list
  done
done
# 2a. Collect counts of all variants impacting any RAS pathway gene
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

  for cancer in PDAC CRAD LUAD ALL; do
    cat $COHORTDIR/data/ras_pathway_coding_vid_lists/$cohort.$cancer.*.nonsyn.vids.list \
    | sort | uniq | wc -l
  done
done | paste - - - - -
# 2b. Union of counts across all genes and cohorts
for cancer in PDAC CRAD LUAD ALL; do
  for cohort in PROFILE HMF TCGA; do
    case $cohort in
      TCGA)
        cat $TCGADIR/data/ras_pathway_coding_vid_lists/$cohort.$cancer.*.nonsyn.vids.list
        ;;
      PROFILE)
        cat $PROFILEDIR/data/ras_pathway_coding_vid_lists/$cohort.$cancer.*.nonsyn.vids.list
        ;;
      HMF)
        cat $HMFDIR/data/ras_pathway_coding_vid_lists/$cohort.$cancer.*.nonsyn.vids.list
        ;;
    esac
  done | sort -V | uniq | wc -l
done | paste <( echo "Union" ) - - - -
# 3a. Compile tables of counts per gene for each cohort & cancer type
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
  while read gene; do
    echo $gene
    for cancer in PDAC CRAD LUAD ALL; do
      cat $COHORTDIR/data/ras_pathway_coding_vid_lists/$cohort.$cancer.$gene.nonsyn.vids.list | wc -l
    done
  done < $CODEDIR/refs/NCI_RAS_pathway.genes.list | paste - - - - - \
  > $COHORTDIR/data/ras_pathway_coding_vid_lists/$cohort.ALL.nonsyn.counts.tsv
done
# 3b. Compile union of counts per gene across cohorts
while read gene; do
  echo $gene
  for cancer in PDAC CRAD LUAD ALL; do
    for cohort in PROFILE HMF TCGA; do
      case $cohort in
        TCGA)
          cat $TCGADIR/data/ras_pathway_coding_vid_lists/$cohort.$cancer.$gene.nonsyn.vids.list
          ;;
        PROFILE)
          cat $PROFILEDIR/data/ras_pathway_coding_vid_lists/$cohort.$cancer.$gene.nonsyn.vids.list
          ;;
        HMF)
          cat $HMFDIR/data/ras_pathway_coding_vid_lists/$cohort.$cancer.$gene.nonsyn.vids.list
          ;;
      esac
    done | sort -V | uniq | wc -l
  done
done < $CODEDIR/refs/NCI_RAS_pathway.genes.list | paste - - - - - \
> $WRKDIR/data/all_cohorts.all_cancers.ras_pathway_genes.nonsyn.counts.tsv
# 4. Format counts of nonsyn variants for top three genes
while read gene; do
  echo -e "\n\n$gene:"
  for cohort in PROFILE HMF TCGA; do
    case $cohort in
      TCGA)
        fgrep -w $gene \
          $TCGADIR/data/ras_pathway_coding_vid_lists/$cohort.ALL.nonsyn.counts.tsv \
        | cut -f2-
        ;;
      PROFILE)
        fgrep -w $gene \
          $PROFILEDIR/data/ras_pathway_coding_vid_lists/$cohort.ALL.nonsyn.counts.tsv \
        | cut -f2-
        ;;
      HMF)
        fgrep -w $gene \
          $HMFDIR/data/ras_pathway_coding_vid_lists/$cohort.ALL.nonsyn.counts.tsv \
        | cut -f2-
        ;;
    esac | paste <( echo "$cohort" ) -
  done
  fgrep -w $gene \
    $WRKDIR/data/all_cohorts.all_cancers.ras_pathway_genes.nonsyn.counts.tsv \
  | cut -f2- | paste <( echo "Union" ) -
done < <( sort -nrk5,5 $WRKDIR/data/all_cohorts.all_cancers.ras_pathway_genes.nonsyn.counts.tsv \
          | head -n3 | cut -f1 )


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

