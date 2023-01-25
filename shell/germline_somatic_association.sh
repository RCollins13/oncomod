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
export bonf_sig=0.0000003218497
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/variant_set_freqs/filtered data/germline_vcfs \
              data/variant_sets/test_sets results results/assoc_stats \
              results/assoc_stats/single results/assoc_stats/merged \
              results/assoc_stats/merged/filtered results/assoc_stats/meta \
              plots/germline_somatic_assoc plots/germline_somatic_assoc/qq; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Determine LD-independent number of individual germline variants to test per cancer type
module load plink/1.90b3
for cancer in PDAC CRAD LUAD SKCM; do
  # Subset VCFs per cohort to non-rare (AC≥10) variants
  for cohort in TCGA PROFILE; do
    case $cohort in
      TCGA)
        COHORTDIR=$TCGADIR
        sample_name="donors"
        ;;
      PROFILE)
        COHORTDIR=$PROFILEDIR
        sample_name="samples"
        ;;
    esac
    bcftools view \
      --min-ac 10 \
      --samples-file $COHORTDIR/data/sample_info/$cohort.$cancer.$sample_name.list \
      -o $COHORTDIR/data/$cohort.RAS_loci.$cancer.ac10plus.vcf.gz \
      -O z \
      $COHORTDIR/data/$cohort.RAS_loci.vcf.gz
    tabix -p vcf -f $COHORTDIR/data/$cohort.RAS_loci.$cancer.ac10plus.vcf.gz
  done

  # Merge AC≥10 VCFs across cohorts
  bcftools merge \
    -m none \
    -o $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.ac10plus.vcf.gz \
    -O z \
    $TCGADIR/data/TCGA.RAS_loci.$cancer.ac10plus.vcf.gz \
    $PROFILEDIR/data/PROFILE.RAS_loci.$cancer.ac10plus.vcf.gz
  tabix -p vcf -f $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.ac10plus.vcf.gz

  # LD prune with PLINK
  plink \
    --vcf $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.ac10plus.vcf.gz \
    --indep-pairwise 100kb 5 0.2 \
    --recode vcf bgz \
    --out $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.ac10plus.pruned
done

# Once merged & LD-pruned, determine number of variants retained per cancer type per gene 
for cancer in PDAC CRAD LUAD SKCM; do
  while read chrom start end gene; do
    bcftools query \
      -f '%ID\n' -r "$chrom" \
      $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.ac10plus.vcf.gz \
    | fgrep -wf - \
      $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.ac10plus.pruned.prune.in \
    | wc -l
  done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz ) | paste -s - \
  | awk -v OFS="\t" '{ print $0, $1+$2+$3 }'
done | paste -s -


### Filter germline variant sets for RAS genes to determine which have sufficient 
### data to be tested (AC≥10)
for cohort in TCGA PROFILE; do
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.germline.burden_sets.freq.tsv.gz
  zcat $WRKDIR/data/variant_sets/$cohort.germline.burden_sets.tsv.gz \
  | awk '{ if ($4 ~ /,/) print $1 }' \
  | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
  | grep -e '^set_id\|^NRAS_\|^HRAS_\|^KRAS_' \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-ac 10 \
    --min-freq 0 \
    --report-ac \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.germline.burden_sets.freq.ac10plus.tsv.gz
done

# Summarize filtered sets
$CODEDIR/scripts/germline_somatic_assoc/summarize_germline_burden_sets.py \
  --burden-sets $WRKDIR/data/variant_set_freqs/filtered/TCGA.germline.burden_sets.freq.ac10plus.tsv.gz \
  --burden-sets $WRKDIR/data/variant_set_freqs/filtered/PROFILE.germline.burden_sets.freq.ac10plus.tsv.gz \
  --out-prefix $WRKDIR/data/variant_sets/test_sets/

# Supplement filtered sets with individual variant IDs per gene & cancer type
for cancer in PDAC CRAD LUAD SKCM; do
  while read chrom start end gene; do
    bcftools query \
      --format '%ID\n' \
      --regions $chrom \
      $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.ac10plus.vcf.gz \
    >> $WRKDIR/data/variant_sets/test_sets/$cancer.$gene.germline_sets.tsv
  done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
done


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
#    (Require mutation sets to have two or more individual mutations, otherwise 
#     they would already be captured by single variant tests)
for cohort in TCGA PROFILE; do
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.somatic.burden_sets.freq.tsv.gz
  zcat $WRKDIR/data/variant_sets/$cohort.somatic.burden_sets.tsv.gz \
  | awk '{ if ($4 ~ /,/) print $1 }' \
  | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
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
for cohort in TCGA PROFILE; do
  $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv $WRKDIR/data/variant_set_freqs/$cohort.somatic.ras_plus_nonRas_comutations.freq.tsv.gz \
    --min-freq 0.05 \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.ras_plus_nonRas_comutations.freq.5pct.tsv.gz
done


### Summarize somatic conditions to test as endpoints for association
$CODEDIR/scripts/germline_somatic_assoc/summarize_somatic_endpoints.py \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.coding_variants.freq.1pct.tsv.gz \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.coding_variants.freq.1pct.tsv.gz \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.other_variants.freq.1pct.tsv.gz \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.other_variants.freq.1pct.tsv.gz \
  --codons $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.recurrently_mutated_codons.freq.1pct.tsv.gz \
  --codons $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.recurrently_mutated_codons.freq.1pct.tsv.gz \
  --burden-sets $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.burden_sets.freq.1pct.tsv.gz \
  --burden-sets $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.burden_sets.freq.1pct.tsv.gz \
  --comutations $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.comutations.freq.1pct.tsv.gz \
  --comutations $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.comutations.freq.1pct.tsv.gz \
  --ras-nonras-comut $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.ras_plus_nonRas_comutations.freq.5pct.tsv.gz \
  --ras-nonras-comut $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.ras_plus_nonRas_comutations.freq.5pct.tsv.gz \
  --transcript-info $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
  --out-prefix $WRKDIR/data/variant_sets/test_sets/


### Annotate somatic and germline endpoint/test sets with their constitutent variant IDs
for context in germline somatic; do
  case $context in
    "germline")
      suffix="sets"
      ;;
    "somatic")
      suffix="endpoints"
      ;;
  esac
  for cancer in PDAC CRAD LUAD SKCM; do
    while read chrom start end gene; do
      $CODEDIR/scripts/data_processing/add_variant_set_members.py \
        --set-list $WRKDIR/data/variant_sets/test_sets/$cancer.$gene.${context}_$suffix.tsv \
        --memberships $WRKDIR/data/variant_sets/PROFILE.$context.burden_sets.tsv.gz \
        --memberships $WRKDIR/data/variant_sets/PROFILE.$context.collapsed_coding_csqs.tsv.gz \
        --memberships $WRKDIR/data/variant_sets/PROFILE.$context.other_single_variants.tsv.gz \
        --memberships $WRKDIR/data/variant_sets/PROFILE.$context.recurrently_mutated_codons.tsv.gz \
        --memberships $WRKDIR/data/variant_sets/TCGA.$context.burden_sets.tsv.gz \
        --memberships $WRKDIR/data/variant_sets/TCGA.$context.collapsed_coding_csqs.tsv.gz \
        --memberships $WRKDIR/data/variant_sets/TCGA.$context.other_single_variants.tsv.gz \
        --memberships $WRKDIR/data/variant_sets/TCGA.$context.recurrently_mutated_codons.tsv.gz
    done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
  done
done


### Submit germline-somatic association jobs
# Ensure most recent version of RASMod R package is installed from source
Rscript -e "install.packages('$CODEDIR/src/RASMod_0.1.tar.gz', \
                             lib='~/R/x86_64-pc-linux-gnu-library/3.6', \
                             type='source', repos=NULL)"
# One submission per (gene, cancer type, cohort)
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      mem=16000
      cpu=4
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      mem=32000
      cpu=8
      ;;
  esac
  for cancer in PDAC CRAD LUAD SKCM; do
    while read chrom start end gene; do
      if ! [ -e $WRKDIR/results/assoc_stats/single/$cohort.$cancer.$gene.sumstats.tsv.gz ] || \
         ! [ -s $WRKDIR/results/assoc_stats/single/$cohort.$cancer.$gene.sumstats.tsv.gz ]; then
        cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_assoc_${cohort}_${cancer}_${gene}.sh
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.single.R \
  --sample-metadata $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
  --cancer-type $cancer \
  --somatic-ad $COHORTDIR/data/$cohort.somatic_variants.dosage.tsv.gz \
  --germline-ad $COHORTDIR/data/$cohort.RAS_loci.dosage.tsv.gz \
  --somatic-variant-sets $WRKDIR/data/variant_sets/test_sets/$cancer.$gene.somatic_endpoints.tsv \
  --germline-variant-sets $WRKDIR/data/variant_sets/test_sets/$cancer.$gene.germline_sets.tsv \
  --outfile $WRKDIR/results/assoc_stats/single/$cohort.$cancer.$gene.sumstats.tsv
gzip -f $WRKDIR/results/assoc_stats/single/$cohort.$cancer.$gene.sumstats.tsv
EOF
        chmod a+x $WRKDIR/LSF/scripts/germline_somatic_assoc_${cohort}_${cancer}_${gene}.sh
        for suf in err log; do
          logfile=$WRKDIR/LSF/logs/germline_somatic_assoc_${cohort}_${cancer}_${gene}.$suf
          if [ -e $logfile ]; then rm $logfile; fi
        done
        bsub -q big-multi -sla miket_sc -R "rusage[mem=${mem}]" -n $cpu \
          -J germline_somatic_assoc_${cohort}_${cancer}_${gene} \
          -o $WRKDIR/LSF/logs/germline_somatic_assoc_${cohort}_${cancer}_${gene}.log \
          -e $WRKDIR/LSF/logs/germline_somatic_assoc_${cohort}_${cancer}_${gene}.err \
          $WRKDIR/LSF/scripts/germline_somatic_assoc_${cohort}_${cancer}_${gene}.sh
      fi
    done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
  done
done


### Plot Q-Qs for each set of association statistics
# 1. Ensure newest version of rCNV2 R package is loaded
cd $WRKDIR/../code/rCNV2 && \
git pull && \
cd - && \
Rscript -e "install.packages('$WRKDIR/../code/rCNV2/source/rCNV2_1.0.1.tar.gz', \
                             lib='~/R/x86_64-pc-linux-gnu-library/3.6', \
                             type='source', repos=NULL)"
# 2. Combine association statistics across all genes per cancer type & cohort
for cohort in TCGA PROFILE; do
  for cancer in PDAC CRAD LUAD SKCM; do
    if [ $( find $WRKDIR/results/assoc_stats/single/ -name "$cohort.$cancer.*.sumstats.tsv.gz" | wc -l ) -eq 3 ]; then
      zcat $WRKDIR/results/assoc_stats/single/$cohort.$cancer.*.sumstats.tsv.gz \
      | grep -v "^#" | sort -nk9,9 \
      | cat <( zcat $WRKDIR/results/assoc_stats/single/$cohort.$cancer.NRAS.sumstats.tsv.gz | head -n1 ) - \
      | gzip -c \
      > $WRKDIR/results/assoc_stats/merged/$cohort.$cancer.sumstats.tsv.gz
    fi
  done
done
# 3. For QQ visualization, restrict to qualifying events in each cohort (based on somatic frequency)
for cohort in TCGA PROFILE; do
  for cancer in PDAC CRAD LUAD SKCM; do
    stats=$WRKDIR/results/assoc_stats/merged/$cohort.$cancer.sumstats.tsv.gz
    if [ -s $stats ]; then
      cat \
        <( zcat $stats | head -n1 ) \
        <( zcat $stats | grep -ve '^COMUT\|^#' | awk '{ if ($4/$3 >= 0.01) print }' ) \
        <( zcat $stats \
           | fgrep -wf \
               <( zcat $WRKDIR/data/variant_set_freqs/$cohort.somatic.gene_comutations.freq.tsv.gz \
                  | sed '1d' | cut -f1 ) \
           | awk '{ if ($4/$3 >= 0.01) print }' ) \
        <( zcat $stats \
           | fgrep -wf \
               <( zcat $WRKDIR/data/variant_set_freqs/$cohort.somatic.ras_plus_nonRas_comutations.freq.tsv.gz \
                  | sed '1d' | cut -f1 ) \
           | awk '{ if ($4/$3 >= 0.05) print }' ) \
      | gzip -c \
      > $WRKDIR/results/assoc_stats/merged/filtered/$cohort.$cancer.sumstats.filtered.tsv.gz
    fi
  done
done
# 4. Plot one QQ for each cancer type & cohort
for cohort in TCGA PROFILE; do
  case $cohort in
    PROFILE)
      alt_cohort=DFCI
      ;;
    *)
      alt_cohort=$cohort
      ;;
  esac
  for cancer in PDAC CRAD LUAD SKCM; do
    stats=$WRKDIR/results/assoc_stats/merged/filtered/$cohort.$cancer.sumstats.filtered.tsv.gz
    if [ -s $stats ]; then
      bsub -q short -sla miket_sc -J plot_qq_single_${cohort}_${cancer} \
        -o $WRKDIR/LSF/logs/plot_qq_single_${cohort}_${cancer}.log \
        -e $WRKDIR/LSF/logs/plot_qq_single_${cohort}_${cancer}.err \
        "$CODEDIR/utils/plot_qq.R \
           --stats $stats \
           --outfile $WRKDIR/plots/germline_somatic_assoc/qq/$cohort.$cancer.qq.png \
           --cancer $cancer \
           --cohort $alt_cohort \
           --p-threshold $bonf_sig"
    fi
  done
done


### Submit meta-analyses for overlapping variants between cohorts
# One submission per cancer type
for cancer in PDAC CRAD LUAD SKCM; do
  TCGA_stats=$WRKDIR/results/assoc_stats/merged/filtered/TCGA.$cancer.sumstats.filtered.tsv.gz
  PROFILE_stats=$WRKDIR/results/assoc_stats/merged/filtered/PROFILE.$cancer.sumstats.filtered.tsv.gz
  if [ -s $TCGA_stats ] && \
     [ -s $PROFILE_stats ]; then
    cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_meta_$cancer.sh
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.meta.R \
  --stats $TCGA_stats \
  --name TCGA \
  --stats $PROFILE_stats \
  --name DFCI \
  --outfile $WRKDIR/results/assoc_stats/meta/$cancer.meta.sumstats.tsv
gzip -f $WRKDIR/results/assoc_stats/meta/$cancer.meta.sumstats.tsv
EOF
    chmod a+x $WRKDIR/LSF/scripts/germline_somatic_meta_$cancer.sh
    for suf in err log; do
      logfile=$WRKDIR/LSF/logs/germline_somatic_meta_$cancer.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
      -J germline_somatic_meta_$cancer \
      -o $WRKDIR/LSF/logs/germline_somatic_meta_$cancer.log \
      -e $WRKDIR/LSF/logs/germline_somatic_meta_$cancer.err \
      $WRKDIR/LSF/scripts/germline_somatic_meta_$cancer.sh
  fi
done
# Once complete, plot one QQ for each cancer type
for cancer in PDAC CRAD LUAD SKCM; do
  if [ -e $WRKDIR/results/assoc_stats/meta/$cancer.meta.sumstats.tsv.gz ]; then
    # Plot one QQ of all sumstats
    $CODEDIR/utils/plot_qq.R \
      --stats $WRKDIR/results/assoc_stats/meta/$cancer.meta.sumstats.tsv.gz \
      --outfile $WRKDIR/plots/germline_somatic_assoc/qq/$cancer.meta_plus_single.qq.png \
      --cancer $cancer \
      --cohort "All Results" \
      --p-threshold $bonf_sig
    # Plot a second QQ of only meta-analyzed sumstats
    zcat $WRKDIR/results/assoc_stats/meta/$cancer.meta.sumstats.tsv.gz \
    | awk -v FS="\t" '{ if ($6>1) print }' | gzip -c \
    > $TMPDIR/$cancer.meta.sumstats.meta_only.tsv.gz
    $CODEDIR/utils/plot_qq.R \
      --stats $TMPDIR/$cancer.meta.sumstats.meta_only.tsv.gz \
      --outfile $WRKDIR/plots/germline_somatic_assoc/qq/$cancer.meta_only.qq.png \
      --cancer $cancer \
      --cohort "Meta-analysis" \
      --p-threshold $bonf_sig
  fi
done


### Meta-analyze results across cancers
if [ -s $WRKDIR/results/assoc_stats/merged/TCGA.PDAC.sumstats.tsv.gz ] && \
   [ -s $WRKDIR/results/assoc_stats/merged/PROFILE.PDAC.sumstats.tsv.gz ] && \
   [ -s $WRKDIR/results/assoc_stats/merged/TCGA.CRAD.sumstats.tsv.gz ] && \
   [ -s $WRKDIR/results/assoc_stats/merged/PROFILE.CRAD.sumstats.tsv.gz ] && \
   [ -s $WRKDIR/results/assoc_stats/merged/TCGA.LUAD.sumstats.tsv.gz ] && \
   [ -s $WRKDIR/results/assoc_stats/merged/PROFILE.LUAD.sumstats.tsv.gz ] && \
   [ -s $WRKDIR/results/assoc_stats/merged/TCGA.SKCM.sumstats.tsv.gz ] && \
   [ -s $WRKDIR/results/assoc_stats/merged/PROFILE.SKCM.sumstats.tsv.gz ]; then
  cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_meta_PanCancer.sh
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.meta.R \
  --stats $WRKDIR/results/assoc_stats/meta/PDAC.meta.sumstats.tsv.gz \
  --name PDAC \
  --stats $WRKDIR/results/assoc_stats/meta/CRAD.meta.sumstats.tsv.gz \
  --name CRAD \
  --stats $WRKDIR/results/assoc_stats/meta/LUAD.meta.sumstats.tsv.gz \
  --name LUAD \
  --stats $WRKDIR/results/assoc_stats/meta/SKCM.meta.sumstats.tsv.gz \
  --name SKCM \
  --outfile $WRKDIR/results/assoc_stats/meta/PanCancer.meta.sumstats.tsv
gzip -f $WRKDIR/results/assoc_stats/meta/PanCancer.meta.sumstats.tsv
EOF
  chmod a+x $WRKDIR/LSF/scripts/germline_somatic_meta_PanCancer.sh
  for suf in err log; do
    logfile=$WRKDIR/LSF/logs/germline_somatic_meta_PanCancer.$suf
    if [ -e $logfile ]; then rm $logfile; fi
  done
  bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
    -J germline_somatic_meta_PanCancer \
    -o $WRKDIR/LSF/logs/germline_somatic_meta_PanCancer.log \
    -e $WRKDIR/LSF/logs/germline_somatic_meta_PanCancer.err \
    $WRKDIR/LSF/scripts/germline_somatic_meta_PanCancer.sh
fi
# Once complete, plot QQ
for cancer in PDAC CRAD LUAD SKCM; do
  if [ -e $WRKDIR/results/assoc_stats/meta/PanCancer.meta.sumstats.tsv.gz ]; then
    $CODEDIR/utils/plot_qq.R \
      --stats $WRKDIR/results/assoc_stats/meta/PanCancer.meta.sumstats.tsv.gz \
      --outfile $WRKDIR/plots/germline_somatic_assoc/qq/PanCancer.meta.qq.png \
      --p-threshold $bonf_sig
  fi
done
