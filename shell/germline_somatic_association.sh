#!/usr/bin/env bash

################################
#    EGFR Modifiers Project    #
################################

# Copyright (c) 2023-Present Ryan L. Collins, and the Gusev/Van Allen Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Germline-somatic interaction analyses for EGFR in LUAD

# Note: intended to be executed on the MGB ERISOne cluster


##################
### Basic prep ###
##################

### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export WRKDIR=/data/gusev/USERS/rlc47/EGFR_modifier_analysis
export CODEDIR=$WRKDIR/../EGFR_modifier_analysis/code/ras_modifiers
export bonf_sig=0.000006924249
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
    --samples-file $COHORTDIR/data/sample_info/$cohort.LUAD.$sample_name.list \
    $COHORTDIR/data/$cohort.EGFR_loci.vcf.gz \
  | bcftools +fill-tags - -- -t AC,AN,F_MISSING,HWE \
  | bcftools view \
    --include 'AC >= 10 & (AN-AC) >= 10 & F_MISSING < 0.5 & HWE>0.000001' \
    -o $COHORTDIR/data/$cohort.EGFR_loci.LUAD.qc_pass.vcf.gz \
    -O z
  tabix -p vcf -f $COHORTDIR/data/$cohort.EGFR_loci.LUAD.qc_pass.vcf.gz
done

# Merge QC-pass VCFs across cohorts
bcftools merge \
  -m none \
  -o $WRKDIR/data/germline_vcfs/all_cohorts.EGFR_loci.LUAD.qc_pass.vcf.gz \
  -O z \
  $TCGADIR/data/TCGA.EGFR_loci.LUAD.qc_pass.vcf.gz \
  $PROFILEDIR/data/PROFILE.EGFR_loci.LUAD.qc_pass.vcf.gz
tabix -p vcf -f $WRKDIR/data/germline_vcfs/all_cohorts.EGFR_loci.LUAD.qc_pass.vcf.gz

# LD prune with PLINK
plink \
  --vcf $WRKDIR/data/germline_vcfs/all_cohorts.EGFR_loci.LUAD.qc_pass.vcf.gz \
  --indep-pairwise 100kb 5 0.2 \
  --recode vcf bgz \
  --out $WRKDIR/data/germline_vcfs/all_cohorts.EGFR_loci.LUAD.qc_pass.pruned

# Once merged & LD-pruned, determine number of variants retained per cancer type per gene 
while read chrom start end gene; do
  bcftools query \
    -f '%ID\n' -r "$chrom" \
    $WRKDIR/data/germline_vcfs/all_cohorts.EGFR_loci.LUAD.qc_pass.vcf.gz \
  | fgrep -wf - \
    $WRKDIR/data/germline_vcfs/all_cohorts.EGFR_loci.LUAD.qc_pass.pruned.prune.in \
  | wc -l
done < <( zcat $WRKDIR/../refs/EGFR_gene.bed.gz ) | paste -s - \
| awk -v OFS="\t" '{ print $0, $1+$2+$3 }'


### Filter germline variant sets for EGFR to determine which have sufficient 
### data to be tested (AC≥10)
for cohort in TCGA PROFILE; do
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.germline.burden_sets.freq.tsv.gz
  zcat $WRKDIR/data/variant_sets/$cohort.germline.burden_sets.tsv.gz \
  | awk '{ if ($4 ~ /,/) print $1 }' \
  | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
  | grep -e '^set_id\|^EGFR_' \
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
while read chrom start end gene; do
  bcftools query \
    --format '%ID\n' \
    --regions $chrom \
    $WRKDIR/data/germline_vcfs/all_cohorts.EGFR_loci.LUAD.qc_pass.vcf.gz \
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
  >> $WRKDIR/data/variant_sets/test_sets/LUAD.$gene.germline_sets.tsv
done < <( zcat $WRKDIR/../refs/EGFR_gene.bed.gz )

# Make lists of PRS to test per cancer type in PROFILE
awk -v cancer=LUAD -v OFS="\t" '{ if ($1==cancer) print $2, $2 }' \
  $CODEDIR/refs/PROFILE_selected_PRS.tsv \
| sort -V \
> $WRKDIR/data/variant_sets/test_sets/LUAD.germline_PRS.tsv


### Filter somatic variants to define list of conditions to test
# 1. Frequent EGFR mutations
for cohort in TCGA PROFILE; do
  # Coding (collapsed by consequence)
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.somatic.coding_variants.freq.tsv.gz
  zcat $WRKDIR/data/variant_sets/$cohort.somatic.collapsed_coding_csqs.tsv.gz \
  | fgrep 'EGFR' | cut -f1 \
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
    -R $WRKDIR/../refs/EGFR_gene.bed.gz \
    $COHORTDIR/data/$cohort.somatic_variants.anno.clean.vcf.gz \
  | cut -f3 | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-freq 0.01 \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.other_variants.freq.1pct.tsv.gz
done

# 2. Recurrently mutated EGFR codons
for cohort in TCGA PROFILE; do
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.somatic.recurrently_mutated_codons.freq.tsv.gz
  zcat $WRKDIR/data/variant_sets/$cohort.somatic.recurrently_mutated_codons.tsv.gz \
  | fgrep -wf \
    <( zcat $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
       | awk -v FS="\t" '{ if ($3 == "EGFR") print $1 }' \
       | cut -f1 -d\. ) \
  | cut -f1 | fgrep -wf - <( zcat $freqs ) | cat <( zcat $freqs | head -n1 ) - \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-freq 0.01 \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.recurrently_mutated_codons.freq.1pct.tsv.gz
done

# 3. Recurrently mutated EGFR exons
for cohort in TCGA PROFILE; do
  freqs=$WRKDIR/data/variant_set_freqs/$cohort.somatic.recurrently_mutated_exons.freq.tsv.gz
  zcat $WRKDIR/data/variant_sets/$cohort.somatic.recurrently_mutated_exons.tsv.gz \
  | fgrep -wf \
    <( zcat $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
       | awk -v FS="\t" '{ if ($3 == "EGFR") print $1 }' \
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
  | grep -e '^set_id\|^EGFR_' \
  | $CODEDIR/scripts/data_processing/filter_freq_table.py \
    --freq-tsv stdin \
    --min-freq 0.01 \
    --outfile $WRKDIR/data/variant_set_freqs/filtered/$cohort.somatic.burden_sets.freq.1pct.tsv.gz
done

# 5. Frequent intra-gene co-mutation pairs involving EGFR
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
    -R $WRKDIR/../refs/EGFR_gene.bed.gz \
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
while read chrom start end gene; do
  $CODEDIR/../GenomicsToolbox/evenSplitter.R \
    -L 750 \
    --shuffle \
    $WRKDIR/data/variant_sets/test_sets/LUAD.$gene.germline_sets.tsv \
    $WRKDIR/data/variant_sets/test_sets/shards/LUAD.$gene.germline_sets.shard_
done < <( zcat $WRKDIR/../refs/EGFR_gene.bed.gz )


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
  while read chrom start end gene; do
    cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_assoc_${cohort}_${cancer}_${gene}.sharded.sh
#!/usr/bin/env bash
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.single.R \
  --sample-metadata $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
  --cancer-type LUAD \
  --somatic-ad $COHORTDIR/data/$cohort.somatic_variants.dosage.tsv.gz \
  --germline-ad $COHORTDIR/data/$cohort.EGFR_loci.dosage.tsv.gz \
  --somatic-variant-sets $WRKDIR/data/variant_sets/test_sets/LUAD.$gene.somatic_endpoints.tsv \
  --germline-variant-sets $WRKDIR/data/variant_sets/test_sets/shards/LUAD.$gene.germline_sets.shard_\$1 \
  --eligible-controls $COHORTDIR/data/sample_info/$cohort.ALL.eligible_controls.list \
  --outfile $WRKDIR/results/assoc_stats/single/$cohort.LUAD.$gene.sumstats.\$1.tsv
gzip -f $WRKDIR/results/assoc_stats/single/$cohort.LUAD.$gene.sumstats.\$1.tsv
EOF
    chmod a+x $WRKDIR/LSF/scripts/germline_somatic_assoc_${cohort}_${cancer}_${gene}.sharded.sh
    n_shards=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                  -name "LUAD.$gene.germline_sets.shard_*" | wc -l )
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
  done < <( zcat $WRKDIR/../refs/EGFR_gene.bed.gz )
done
# Submit one PRS association job per cancer type per gene (for PROFILE only)
prs_sets=$WRKDIR/data/variant_sets/test_sets/LUAD.germline_PRS.tsv
if [ $( cat $prs_sets | wc -l ) -gt 0 ]; then
  while read chrom start end gene; do
    cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_assoc_PROFILE_${cancer}_${gene}.PRS.sh
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.single.R \
  --sample-metadata $PROFILEDIR/data/sample_info/PROFILE.ALL.sample_metadata.tsv.gz \
  --cancer-type LUAD \
  --somatic-ad $PROFILEDIR/data/PROFILE.somatic_variants.dosage.tsv.gz \
  --germline-ad $PROFILEDIR/data/PROFILE.PRS.tsv.gz \
  --somatic-variant-sets $WRKDIR/data/variant_sets/test_sets/LUAD.$gene.somatic_endpoints.tsv \
  --germline-variant-sets $prs_sets \
  --normalize-germline-ad \
  --eligible-controls $PROFILEDIR/data/sample_info/PROFILE.ALL.eligible_controls.list \
  --outfile $WRKDIR/results/assoc_stats/single/PROFILE.LUAD.$gene.sumstats.PRS.tsv
gzip -f $WRKDIR/results/assoc_stats/single/PROFILE.LUAD.$gene.sumstats.PRS.tsv
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
  done < <( zcat $WRKDIR/../refs/EGFR_gene.bed.gz )
fi
# Find missing/incomplete shards
for cohort in TCGA PROFILE; do
  while read chrom start end gene; do
    n_shards=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                    -name "LUAD.$gene.germline_sets.shard_*" | wc -l )
    for i in $( seq 1 $n_shards ); do
      if ! [ -s $WRKDIR/results/assoc_stats/single/$cohort.LUAD.$gene.sumstats.$i.tsv.gz ]; then
        echo -e "$cohort\tLUAD\t$gene\t$i"
      fi
    done 
    if [ $cohort == "PROFILE" ] && \
       [ -s $WRKDIR/data/variant_sets/test_sets/LUAD.germline_PRS.tsv ]; then
      if ! [ -s $WRKDIR/results/assoc_stats/single/PROFILE.LUAD.$gene.sumstats.PRS.tsv.gz ]; then
        echo -e "$cohort\tLUAD\t$gene\tPRS"
      fi
    fi
  done < <( zcat $WRKDIR/../refs/EGFR_gene.bed.gz )
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
  n_expected=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                  -name "LUAD.*.germline_sets.shard_*" | wc -l )
  if [ $cohort == "PROFILE" ] && \
     [ -s $WRKDIR/data/variant_sets/test_sets/LUAD.germline_PRS.tsv ]; then
    n_prs_expected=3
  else
    n_prs_expected=0
  fi
  n_complete=$( find $WRKDIR/results/assoc_stats/single/ \
                  -name "$cohort.LUAD.*.sumstats.*.tsv.gz" | wc -l )
  if [ $(( $n_expected + $n_prs_expected )) -eq $n_complete ]; then
    zcat $WRKDIR/results/assoc_stats/single/$cohort.LUAD.*.sumstats.*.tsv.gz \
    | grep -v "^#" | sort -nrk9,9 \
    | cat <( zcat $WRKDIR/results/assoc_stats/single/$cohort.LUAD.EGFR.sumstats.1.tsv.gz | head -n1 ) - \
    | gzip -c \
    > $WRKDIR/results/assoc_stats/merged/$cohort.LUAD.sumstats.tsv.gz
  fi
done
# 3. For QQ visualization, restrict to qualifying events in each cohort (based on somatic frequency)
for cohort in TCGA PROFILE; do
  stats=$WRKDIR/results/assoc_stats/merged/$cohort.LUAD.sumstats.tsv.gz
  if [ -s $stats ]; then
    cat \
      <( zcat $stats | head -n1 ) \
      <( zcat $stats | grep -ve '^#' | awk '{ if ($4/$3 >= 0.01) print }' ) \
    | gzip -c \
    > $WRKDIR/results/assoc_stats/merged/filtered/$cohort.LUAD.sumstats.filtered.tsv.gz
  fi
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
  stats=$WRKDIR/results/assoc_stats/merged/filtered/$cohort.LUAD.sumstats.filtered.tsv.gz
  if [ -s $stats ]; then
    bsub -q short -sla miket_sc -J plot_qq_single_${cohort}_${cancer} \
      -o $WRKDIR/LSF/logs/plot_qq_single_${cohort}_${cancer}.log \
      -e $WRKDIR/LSF/logs/plot_qq_single_${cohort}_${cancer}.err \
      "$CODEDIR/utils/plot_qq.R \
         --stats $stats \
         --outfile $WRKDIR/plots/germline_somatic_assoc/qq/$cohort.LUAD.qq.png \
         --cancer LUAD \
         --cohort $alt_cohort \
         --p-threshold $bonf_sig"
  fi
done


### Mega-analysis (pooled) of all cohorts per cancer type
# One submission per gene & cancer type
while read chrom start end gene; do
  cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_${cancer}_${gene}.sharded.sh
#!/usr/bin/env bash
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.pooled.R \
--sample-metadata $TCGADIR/data/sample_info/TCGA.ALL.sample_metadata.tsv.gz \
--somatic-ad $TCGADIR/data/TCGA.somatic_variants.dosage.tsv.gz \
--germline-ad $TCGADIR/data/TCGA.EGFR_loci.dosage.tsv.gz \
--name TCGA \
--eligible-controls $TCGADIR/data/sample_info/TCGA.ALL.eligible_controls.list \
--sample-metadata $PROFILEDIR/data/sample_info/PROFILE.ALL.sample_metadata.tsv.gz \
--somatic-ad $PROFILEDIR/data/PROFILE.somatic_variants.dosage.tsv.gz \
--germline-ad $PROFILEDIR/data/PROFILE.EGFR_loci.dosage.tsv.gz \
--name DFCI \
--eligible-controls $PROFILEDIR/data/sample_info/PROFILE.ALL.eligible_controls.list \
--cancer-type LUAD \
--somatic-variant-sets $WRKDIR/data/variant_sets/test_sets/LUAD.$gene.somatic_endpoints.tsv \
--germline-variant-sets $WRKDIR/data/variant_sets/test_sets/shards/LUAD.$gene.germline_sets.shard_\$1 \
--outfile $WRKDIR/results/assoc_stats/single/pooled.LUAD.$gene.sumstats.\$1.tsv
gzip -f $WRKDIR/results/assoc_stats/single/pooled.LUAD.$gene.sumstats.\$1.tsv
EOF
  chmod a+x $WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_${cancer}_${gene}.sharded.sh
  n_shards=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                -name "LUAD.$gene.germline_sets.shard_*" | wc -l )
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
done < <( zcat $WRKDIR/../refs/EGFR_gene.bed.gz )
# Find missing shards
while read chrom start end gene; do
  n_shards=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                  -name "LUAD.$gene.germline_sets.shard_*" | wc -l )
  for i in $( seq 1 $n_shards ); do
    if ! [ -s $WRKDIR/results/assoc_stats/single/pooled.LUAD.$gene.sumstats.$i.tsv.gz ]; then
      echo -e "LUAD\t$gene\t$i"
    fi
  done 
done < <( zcat $WRKDIR/../refs/EGFR_gene.bed.gz )
# Once complete, pool results per cancer type
n_expected=$( find $WRKDIR/data/variant_sets/test_sets/shards/ \
                -name "LUAD.*.germline_sets.shard_*" | wc -l )
n_complete=$( find $WRKDIR/results/assoc_stats/single/ \
                -name "pooled.LUAD.*.sumstats.*.tsv.gz" | wc -l )
if [ $n_expected -eq $n_complete ]; then
  zcat \
    $WRKDIR/results/assoc_stats/single/pooled.LUAD.*.sumstats.*.tsv.gz \
    $WRKDIR/results/assoc_stats/single/PROFILE.LUAD.*.sumstats.PRS.tsv.gz \
  | grep -v "^#" | sort -nrk9,9 \
  | cat <( zcat $WRKDIR/results/assoc_stats/single/pooled.LUAD.EGFR.sumstats.1.tsv.gz | head -n1 ) - \
  | gzip -c \
  > $WRKDIR/results/assoc_stats/merged/pooled.LUAD.sumstats.tsv.gz
fi
# Once complete, plot one QQ for each cancer type
stats=$WRKDIR/results/assoc_stats/merged/pooled.LUAD.sumstats.tsv.gz
if [ -s $stats ]; then
  bsub -q short -sla miket_sc -J plot_qq_pooled_LUAD \
    -o $WRKDIR/LSF/logs/plot_qq_pooled_LUAD.log \
    -e $WRKDIR/LSF/logs/plot_qq_pooled_LUAD.err \
    "$CODEDIR/utils/plot_qq.R \
       --stats $stats \
       --outfile $WRKDIR/plots/germline_somatic_assoc/qq/LUAD.pooled.qq.png \
       --cancer LUAD \
       --cohort \"Pooled Analysis\" \
       --p-threshold $bonf_sig"
fi


### Gather significant hits in any individual cohort or pooled analysis
$CODEDIR/scripts/germline_somatic_assoc/get_sig_hits.py \
  --sumstats $WRKDIR/results/assoc_stats/merged/pooled.LUAD.sumstats.tsv.gz \
  --cohort-name Pooled \
  --sumstats $WRKDIR/results/assoc_stats/merged/PROFILE.LUAD.sumstats.tsv.gz \
  --cohort-name PROFILE \
  --sumstats $WRKDIR/results/assoc_stats/merged/TCGA.LUAD.sumstats.tsv.gz \
  --cohort-name TCGA \
  --p-cutoff $bonf_sig

