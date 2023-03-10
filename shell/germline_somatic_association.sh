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
export bonf_sig=0.0000006155892
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/variant_set_freqs/filtered data/germline_vcfs \
              data/variant_sets/test_sets data/variant_sets/test_sets/shards \
              results results/assoc_stats results/assoc_stats/single \
              results/assoc_stats/merged results/assoc_stats/merged/filtered \
              results/assoc_stats/meta plots/germline_somatic_assoc \
              plots/germline_somatic_assoc/qq; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Determine LD-independent number of individual germline variants to test per cancer type
module load plink/1.90b3
for cancer in PDAC CRAD LUAD SKCM; do
  # Subset VCFs per cohort to:
  # 1. Non-rare (10≤AC≤(2*N_samples - 10)) variants
  # 2. called in at least 75% of samples and 
  # 3. in HWE
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
      --samples-file $COHORTDIR/data/sample_info/$cohort.$cancer.$sample_name.list \
      $COHORTDIR/data/$cohort.RAS_loci.vcf.gz \
    | bcftools +fill-tags - -- -t AC,AN,F_MISSING,HWE \
    | bcftools view \
      --include 'AC >= 10 & (AN-AC) >= 10 & F_MISSING < 0.5 & HWE>0.000001' \
      -o $COHORTDIR/data/$cohort.RAS_loci.$cancer.qc_pass.vcf.gz \
      -O z
    tabix -p vcf -f $COHORTDIR/data/$cohort.RAS_loci.$cancer.qc_pass.vcf.gz
  done

  # Merge QC-pass VCFs across cohorts
  bcftools merge \
    -m none \
    -o $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.qc_pass.vcf.gz \
    -O z \
    $TCGADIR/data/TCGA.RAS_loci.$cancer.qc_pass.vcf.gz \
    $PROFILEDIR/data/PROFILE.RAS_loci.$cancer.qc_pass.vcf.gz
  tabix -p vcf -f $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.qc_pass.vcf.gz

  # LD prune with PLINK
  plink \
    --vcf $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.qc_pass.vcf.gz \
    --indep-pairwise 100kb 5 0.2 \
    --recode vcf bgz \
    --out $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.qc_pass.pruned
done

# Once merged & LD-pruned, determine number of variants retained per cancer type per gene 
for cancer in PDAC CRAD LUAD SKCM; do
  while read chrom start end gene; do
    bcftools query \
      -f '%ID\n' -r "$chrom" \
      $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.qc_pass.vcf.gz \
    | fgrep -wf - \
      $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.qc_pass.pruned.prune.in \
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
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.germline.burden_sets.freq.qc_pass.tsv.gz
done

# Summarize filtered sets
$CODEDIR/scripts/germline_somatic_assoc/summarize_germline_burden_sets.py \
  --burden-sets $WRKDIR/data/variant_set_freqs/filtered/TCGA.germline.burden_sets.freq.qc_pass.tsv.gz \
  --burden-sets $WRKDIR/data/variant_set_freqs/filtered/PROFILE.germline.burden_sets.freq.qc_pass.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/PROFILE.germline.burden_sets.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/PROFILE.germline.collapsed_coding_csqs.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/PROFILE.germline.other_single_variants.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/PROFILE.germline.recurrently_mutated_codons.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/TCGA.germline.burden_sets.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/TCGA.germline.collapsed_coding_csqs.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/TCGA.germline.other_single_variants.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/TCGA.germline.recurrently_mutated_codons.tsv.gz \
  --out-prefix $WRKDIR/data/variant_sets/test_sets/

# Supplement filtered sets with individual variant IDs per gene & cancer type
for cancer in PDAC CRAD LUAD SKCM; do
  while read chrom start end gene; do
    bcftools query \
      --format '%ID\n' \
      --regions $chrom \
      $WRKDIR/data/germline_vcfs/all_cohorts.RAS_loci.$cancer.qc_pass.vcf.gz \
    | $CODEDIR/scripts/data_processing/add_variant_set_members.py \
      --set-list stdin \
      --memberships $WRKDIR/data/variant_sets/PROFILE.germline.burden_sets.tsv.gz \
      --memberships $WRKDIR/data/variant_sets/PROFILE.germline.collapsed_coding_csqs.tsv.gz \
      --memberships $WRKDIR/data/variant_sets/PROFILE.germline.other_single_variants.tsv.gz \
      --memberships $WRKDIR/data/variant_sets/PROFILE.germline.recurrently_mutated_codons.tsv.gz \
      --memberships $WRKDIR/data/variant_sets/TCGA.germline.burden_sets.tsv.gz \
      --memberships $WRKDIR/data/variant_sets/TCGA.germline.collapsed_coding_csqs.tsv.gz \
      --memberships $WRKDIR/data/variant_sets/TCGA.germline.other_single_variants.tsv.gz \
      --memberships $WRKDIR/data/variant_sets/TCGA.germline.recurrently_mutated_codons.tsv.gz \
    >> $WRKDIR/data/variant_sets/test_sets/$cancer.$gene.germline_sets.tsv
  done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
done

# Make lists of PRS to test per cancer type in PROFILE
for cancer in PDAC CRAD LUAD SKCM; do
  awk -v cancer=$cancer -v OFS="\t" '{ if ($1==cancer) print $2, $2 }' \
    $CODEDIR/refs/PROFILE_selected_PRS.tsv \
  | sort -V \
  > $WRKDIR/data/variant_sets/test_sets/$cancer.germline_PRS.tsv
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

# 3. Recurrently mutated RAS exons
for cohort in TCGA PROFILE; do
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.somatic.recurrently_mutated_exons.freq.tsv.gz
  zcat $WRKDIR/data/variant_sets/$cohort.somatic.recurrently_mutated_exons.tsv.gz \
  | fgrep -wf \
    <( zcat $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
       | awk -v FS="\t" '{ if ($3 == "NRAS" || $3 == "HRAS" || $3 == "KRAS") print $1 }' \
       | cut -f1 -d\. ) \
  | cut -f1 | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-freq 0.01 \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.recurrently_mutated_exons.freq.1pct.tsv.gz
done

# 4. Functional mutation sets
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

# 5. Frequent intra-gene co-mutation pairs involving RAS
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


### Summarize somatic conditions to test as endpoints for association
$CODEDIR/scripts/germline_somatic_assoc/summarize_somatic_endpoints.py \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.coding_variants.freq.1pct.tsv.gz \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.coding_variants.freq.1pct.tsv.gz \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.other_variants.freq.1pct.tsv.gz \
  --mutations $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.other_variants.freq.1pct.tsv.gz \
  --codons $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.recurrently_mutated_codons.freq.1pct.tsv.gz \
  --codons $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.recurrently_mutated_codons.freq.1pct.tsv.gz \
  --exons $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.recurrently_mutated_exons.freq.1pct.tsv.gz \
  --exons $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.recurrently_mutated_exons.freq.1pct.tsv.gz \
  --burden-sets $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.burden_sets.freq.1pct.tsv.gz \
  --burden-sets $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.burden_sets.freq.1pct.tsv.gz \
  --comutations $WRKDIR/data/variant_set_freqs/filtered/TCGA.somatic.comutations.freq.1pct.tsv.gz \
  --comutations $WRKDIR/data/variant_set_freqs/filtered/PROFILE.somatic.comutations.freq.1pct.tsv.gz \
  --transcript-info $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/PROFILE.somatic.burden_sets.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/PROFILE.somatic.collapsed_coding_csqs.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/PROFILE.somatic.other_single_variants.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/PROFILE.somatic.recurrently_mutated_codons.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/PROFILE.somatic.recurrently_mutated_exons.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/TCGA.somatic.burden_sets.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/TCGA.somatic.collapsed_coding_csqs.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/TCGA.somatic.other_single_variants.tsv.gz \
  --memberships $WRKDIR/data/variant_sets/TCGA.somatic.recurrently_mutated_exons.tsv.gz \
  --out-prefix $WRKDIR/data/variant_sets/test_sets/


### Shard germline test sets for improved parallelization
for cancer in PDAC CRAD LUAD SKCM; do
  while read chrom start end gene; do
    $CODEDIR/../GenomicsToolbox/evenSplitter.R \
      -L 750 \
      --shuffle \
      $WRKDIR/data/variant_sets/test_sets/$cancer.$gene.germline_sets.tsv \
      $WRKDIR/data/variant_sets/test_sets/shards/$cancer.$gene.germline_sets.shard_
  done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
done


## Submit germline-somatic association jobs
# Ensure most recent version of RASMod R package is installed from source
Rscript -e "install.packages('$CODEDIR/src/RASMod_0.1.tar.gz', \
                             lib='~/R/x86_64-pc-linux-gnu-library/3.6', \
                             type='source', repos=NULL)"
# One submission per (gene, cancer type, cohort)
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      mem=6000
      queue=normal
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      mem=12000
      queue=big
      ;;
  esac
  for cancer in PDAC CRAD LUAD SKCM; do
    while read chrom start end gene; do
      cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_assoc_${cohort}_${cancer}_${gene}.sharded.sh
#!/usr/bin/env bash
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.single.R \
  --sample-metadata $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
  --cancer-type $cancer \
  --somatic-ad $COHORTDIR/data/$cohort.somatic_variants.dosage.tsv.gz \
  --germline-ad $COHORTDIR/data/$cohort.RAS_loci.dosage.tsv.gz \
  --somatic-variant-sets $WRKDIR/data/variant_sets/test_sets/$cancer.$gene.somatic_endpoints.tsv \
  --germline-variant-sets $WRKDIR/data/variant_sets/test_sets/shards/$cancer.$gene.germline_sets.shard_\$1 \
  --eligible-controls $COHORTDIR/data/sample_info/$cohort.ALL.eligible_controls.list \
  --outfile $WRKDIR/results/assoc_stats/single/$cohort.$cancer.$gene.sumstats.\$1.tsv
gzip -f $WRKDIR/results/assoc_stats/single/$cohort.$cancer.$gene.sumstats.\$1.tsv
EOF
      chmod a+x $WRKDIR/LSF/scripts/germline_somatic_assoc_${cohort}_${cancer}_${gene}.sharded.sh
      n_shards=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                    -name "$cancer.$gene.germline_sets.shard_*" | wc -l )
      for i in $( seq 1 $n_shards ); do 
        for suf in err log; do
          logfile=$WRKDIR/LSF/logs/germline_somatic_assoc_${cohort}_${cancer}_${gene}.$i.$suf
          if [ -e $logfile ]; then rm $logfile; fi
        done
        bsub -q $queue -sla miket_sc -R "rusage[mem=$mem]" \
          -J germline_somatic_assoc_${cohort}_${cancer}_${gene}_$i \
          -o $WRKDIR/LSF/logs/germline_somatic_assoc_${cohort}_${cancer}_${gene}.$i.log \
          -e $WRKDIR/LSF/logs/germline_somatic_assoc_${cohort}_${cancer}_${gene}.$i.err \
          "$WRKDIR/LSF/scripts/germline_somatic_assoc_${cohort}_${cancer}_${gene}.sharded.sh $i"
      done
    done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
  done
done
# Submit one PRS association job per cancer type per gene (for PROFILE only)
for cancer in PDAC CRAD LUAD SKCM; do
  prs_sets=$WRKDIR/data/variant_sets/test_sets/$cancer.germline_PRS.tsv
  if [ $( cat $prs_sets | wc -l ) -gt 0 ]; then
    while read chrom start end gene; do
      cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_assoc_PROFILE_${cancer}_${gene}.PRS.sh
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.single.R \
  --sample-metadata $PROFILEDIR/data/sample_info/PROFILE.ALL.sample_metadata.tsv.gz \
  --cancer-type $cancer \
  --somatic-ad $PROFILEDIR/data/PROFILE.somatic_variants.dosage.tsv.gz \
  --germline-ad $PROFILEDIR/data/PROFILE.PRS.tsv.gz \
  --somatic-variant-sets $WRKDIR/data/variant_sets/test_sets/$cancer.$gene.somatic_endpoints.tsv \
  --germline-variant-sets $prs_sets \
  --normalize-germline-ad \
  --eligible-controls $PROFILEDIR/data/sample_info/PROFILE.ALL.eligible_controls.list \
  --outfile $WRKDIR/results/assoc_stats/single/PROFILE.$cancer.$gene.sumstats.PRS.tsv
gzip -f $WRKDIR/results/assoc_stats/single/PROFILE.$cancer.$gene.sumstats.PRS.tsv
EOF
      chmod a+x $WRKDIR/LSF/scripts/germline_somatic_assoc_PROFILE_${cancer}_${gene}.PRS.sh
      for suf in err log; do
        logfile=$WRKDIR/LSF/logs/germline_somatic_assoc_PROFILE_${cancer}_${gene}.PRS.$suf
        if [ -e $logfile ]; then rm $logfile; fi
      done
      bsub -q big -sla miket_sc -R "rusage[mem=12000]" \
        -J germline_somatic_assoc_PROFILE_${cancer}_${gene}_PRS \
        -o $WRKDIR/LSF/logs/germline_somatic_assoc_PROFILE_${cancer}_${gene}.PRS.log \
        -e $WRKDIR/LSF/logs/germline_somatic_assoc_PROFILE_${cancer}_${gene}.PRS.err \
        $WRKDIR/LSF/scripts/germline_somatic_assoc_PROFILE_${cancer}_${gene}.PRS.sh
    done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
  fi
done
# Find missing/incomplete shards
for cohort in TCGA PROFILE; do
  for cancer in PDAC CRAD LUAD SKCM; do
    while read chrom start end gene; do
      n_shards=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                      -name "$cancer.$gene.germline_sets.shard_*" | wc -l )
      for i in $( seq 1 $n_shards ); do
        if ! [ -s $WRKDIR/results/assoc_stats/single/$cohort.$cancer.$gene.sumstats.$i.tsv.gz ]; then
          echo -e "$cohort\t$cancer\t$gene\t$i"
        fi
      done 
      if [ $cohort == "PROFILE" ] && \
         [ -s $WRKDIR/data/variant_sets/test_sets/$cancer.germline_PRS.tsv ]; then
        if ! [ -s $WRKDIR/results/assoc_stats/single/PROFILE.$cancer.$gene.sumstats.PRS.tsv.gz ]; then
          echo -e "$cohort\t$cancer\t$gene\tPRS"
        fi
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
# 2. Combine association statistics across all shards per cancer type & cohort
for cohort in TCGA PROFILE; do
  for cancer in PDAC CRAD LUAD SKCM; do
    n_expected=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                    -name "$cancer.*.germline_sets.shard_*" | wc -l )
    if [ $cohort == "PROFILE" ] && \
       [ -s $WRKDIR/data/variant_sets/test_sets/$cancer.germline_PRS.tsv ]; then
      n_prs_expected=3
    else
      n_prs_expected=0
    fi
    n_complete=$( find $WRKDIR/results/assoc_stats/single/ \
                    -name "$cohort.$cancer.*.sumstats.*.tsv.gz" | wc -l )
    if [ $(( $n_expected + $n_prs_expected )) -eq $n_complete ]; then
      zcat $WRKDIR/results/assoc_stats/single/$cohort.$cancer.*.sumstats.*.tsv.gz \
      | grep -v "^#" | sort -nrk9,9 \
      | cat <( zcat $WRKDIR/results/assoc_stats/single/$cohort.$cancer.NRAS.sumstats.1.tsv.gz | head -n1 ) - \
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
        <( zcat $stats | grep -ve '^#' | awk '{ if ($4/$3 >= 0.01) print }' ) \
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


### Mega-analysis (pooled) of all cohorts per cancer type
# One submission per gene & cancer type
for cancer in PDAC CRAD LUAD SKCM; do
  while read chrom start end gene; do
    cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_${cancer}_${gene}.sharded.sh
#!/usr/bin/env bash
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.pooled.R \
--sample-metadata $TCGADIR/data/sample_info/TCGA.ALL.sample_metadata.tsv.gz \
--somatic-ad $TCGADIR/data/TCGA.somatic_variants.dosage.tsv.gz \
--germline-ad $TCGADIR/data/TCGA.RAS_loci.dosage.tsv.gz \
--name TCGA \
--eligible-controls $TCGADIR/data/sample_info/TCGA.ALL.eligible_controls.list \
--sample-metadata $PROFILEDIR/data/sample_info/PROFILE.ALL.sample_metadata.tsv.gz \
--somatic-ad $PROFILEDIR/data/PROFILE.somatic_variants.dosage.tsv.gz \
--germline-ad $PROFILEDIR/data/PROFILE.RAS_loci.dosage.tsv.gz \
--name DFCI \
--eligible-controls $PROFILEDIR/data/sample_info/PROFILE.ALL.eligible_controls.list \
--cancer-type $cancer \
--somatic-variant-sets $WRKDIR/data/variant_sets/test_sets/$cancer.$gene.somatic_endpoints.tsv \
--germline-variant-sets $WRKDIR/data/variant_sets/test_sets/shards/$cancer.$gene.germline_sets.shard_\$1 \
--outfile $WRKDIR/results/assoc_stats/single/pooled.$cancer.$gene.sumstats.\$1.tsv
gzip -f $WRKDIR/results/assoc_stats/single/pooled.$cancer.$gene.sumstats.\$1.tsv
EOF
    chmod a+x $WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_${cancer}_${gene}.sharded.sh
    n_shards=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                  -name "$cancer.$gene.germline_sets.shard_*" | wc -l )
    for i in $( seq 1 $n_shards ); do 
      for suf in err log; do
        logfile=$WRKDIR/LSF/logs/germline_somatic_assoc_pooled_${cancer}_${gene}.$i.$suf
        if [ -e $logfile ]; then rm $logfile; fi
      done
      bsub -q big -sla miket_sc -R "rusage[mem=32000]" \
        -J germline_somatic_assoc_pooled_${cancer}_${gene}_$i \
        -o $WRKDIR/LSF/logs/germline_somatic_assoc_pooled_${cancer}_${gene}.$i.log \
        -e $WRKDIR/LSF/logs/germline_somatic_assoc_pooled_${cancer}_${gene}.$i.err \
        "$WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_${cancer}_${gene}.sharded.sh $i"
    done
  done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
done
# Find missing shards
for cancer in PDAC CRAD LUAD SKCM; do
  while read chrom start end gene; do
    n_shards=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                    -name "$cancer.$gene.germline_sets.shard_*" | wc -l )
    for i in $( seq 1 $n_shards ); do
      if ! [ -s $WRKDIR/results/assoc_stats/single/pooled.$cancer.$gene.sumstats.$i.tsv.gz ]; then
        echo -e "$cancer\t$gene\t$i"
      fi
    done 
  done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
done
# Once complete, pool results per cancer type
for cancer in PDAC CRAD LUAD SKCM; do
  n_expected=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                  -name "$cancer.*.germline_sets.shard_*" | wc -l )
  n_complete=$( find $WRKDIR/results/assoc_stats/single/ \
                  -name "pooled.$cancer.*.sumstats.*.tsv.gz" | wc -l )
  if [ $n_expected -eq $n_complete ]; then
    zcat \
      $WRKDIR/results/assoc_stats/single/pooled.$cancer.*.sumstats.*.tsv.gz \
      $WRKDIR/results/assoc_stats/single/PROFILE.$cancer.*.sumstats.PRS.tsv.gz \
    | grep -v "^#" | sort -nrk9,9 \
    | cat <( zcat $WRKDIR/results/assoc_stats/single/pooled.$cancer.NRAS.sumstats.1.tsv.gz | head -n1 ) - \
    | gzip -c \
    > $WRKDIR/results/assoc_stats/merged/pooled.$cancer.sumstats.tsv.gz
  fi
done
# Once complete, plot one QQ for each cancer type
for cancer in PDAC CRAD LUAD SKCM; do
  stats=$WRKDIR/results/assoc_stats/merged/pooled.$cancer.sumstats.tsv.gz
  if [ -s $stats ]; then
    bsub -q short -sla miket_sc -J plot_qq_pooled_$cancer \
      -o $WRKDIR/LSF/logs/plot_qq_pooled_$cancer.log \
      -e $WRKDIR/LSF/logs/plot_qq_pooled_$cancer.err \
      "$CODEDIR/utils/plot_qq.R \
         --stats $stats \
         --outfile $WRKDIR/plots/germline_somatic_assoc/qq/$cancer.pooled.qq.png \
         --cancer $cancer \
         --cohort \"Pooled Analysis\" \
         --p-threshold $bonf_sig"
  fi
done


### Meta-analyze results across cancers from mega-analysis (pooled)
any_missing=0
for cancer in PDAC CRAD LUAD SKCM; do
  if ! [ -s $WRKDIR/results/assoc_stats/merged/pooled.$cancer.sumstats.tsv.gz ]; then
    any_missing=1
  fi
done
if [ $any_missing -eq 0 ]; then
  cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_pooled_meta_PanCancer.sh
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.meta.R \
  --stats $WRKDIR/results/assoc_stats/merged/pooled.PDAC.sumstats.tsv.gz \
  --name PDAC \
  --stats $WRKDIR/results/assoc_stats/merged/pooled.CRAD.sumstats.tsv.gz \
  --name CRAD \
  --stats $WRKDIR/results/assoc_stats/merged/pooled.LUAD.sumstats.tsv.gz \
  --name LUAD \
  --stats $WRKDIR/results/assoc_stats/merged/pooled.SKCM.sumstats.tsv.gz \
  --name SKCM \
  --drop-frequencies \
  --outfile $WRKDIR/results/assoc_stats/meta/PanCancer.pooled.meta.sumstats.tsv
gzip -f $WRKDIR/results/assoc_stats/meta/PanCancer.pooled.meta.sumstats.tsv
EOF
  chmod a+x $WRKDIR/LSF/scripts/germline_somatic_pooled_meta_PanCancer.sh
  for suf in err log; do
    logfile=$WRKDIR/LSF/logs/germline_somatic_pooled_meta_PanCancer.$suf
    if [ -e $logfile ]; then rm $logfile; fi
  done
  bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
    -J germline_somatic_pooled_meta_PanCancer \
    -o $WRKDIR/LSF/logs/germline_somatic_pooled_meta_PanCancer.log \
    -e $WRKDIR/LSF/logs/germline_somatic_pooled_meta_PanCancer.err \
    $WRKDIR/LSF/scripts/germline_somatic_pooled_meta_PanCancer.sh
fi
# Once complete, subset to tests in two or more cancers and plot QQ
stats=$WRKDIR/results/assoc_stats/meta/PanCancer.pooled.meta.sumstats.tsv.gz
if [ -s $stats ]; then
  zcat $stats | awk '{ if ($1 ~ "#" || $3>1) print $0 }' | gzip -c \
  > $TMPDIR/PanCancer.pooled.meta.stats.tsv.gz
  $CODEDIR/utils/plot_qq.R \
    --stats $TMPDIR/PanCancer.pooled.meta.stats.tsv.gz \
    --outfile $WRKDIR/plots/germline_somatic_assoc/qq/PanCancer.pooled.meta.qq.png \
    --p-threshold $bonf_sig
fi
