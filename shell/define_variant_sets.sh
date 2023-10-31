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
export HMFDIR=/data/gusev/USERS/rlc47/HMF
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export CODEDIR=$WRKDIR/../code/oncomod
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/variant_sets data/variant_set_freqs; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Generate maps of somatic variants with the same protein consequence
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
  for context in germline somatic; do
    case $context in
      germline)
        subset="RAS_loci"
        ;;
      somatic)
        subset="somatic_variants"
        ;;
    esac
    for suf in err log; do
      logfile=$WRKDIR/LSF/logs/collapse_coding_csqs_${cohort}_$context.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q normal -sla miket_sc \
      -J collapse_coding_csqs_${cohort}_$context \
      -o $WRKDIR/LSF/logs/collapse_coding_csqs_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/collapse_coding_csqs_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/vep2csqTable.py \
        $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz \
        $WRKDIR/data/variant_sets/$cohort.$context.collapsed_coding_csqs.tsv.gz"
  done
done


### Identify recurrently mutated codons from the output of collapsed coding csqs, above
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
  for context in germline somatic; do
    for suf in err log; do
      logfile=$WRKDIR/LSF/logs/find_recurrent_codons_${cohort}_$context.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q short -sla miket_sc \
      -J find_recurrent_codons_${cohort}_$context \
      -o $WRKDIR/LSF/logs/find_recurrent_codons_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/find_recurrent_codons_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/find_recurrently_mutated_codons.py \
        $WRKDIR/data/variant_sets/$cohort.$context.collapsed_coding_csqs.tsv.gz \
       | gzip -c > $WRKDIR/data/variant_sets/$cohort.$context.recurrently_mutated_codons.tsv.gz"
  done
done


## Identify recurrently mutated exons from the output of collapsed codons, above
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
  for context in germline somatic; do
    for suf in err log; do
      logfile=$WRKDIR/LSF/logs/find_recurrent_exons_${cohort}_$context.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q short \
      -J find_recurrent_exons_${cohort}_$context \
      -o $WRKDIR/LSF/logs/find_recurrent_exons_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/find_recurrent_exons_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/find_recurrently_mutated_exons.py \
        $WRKDIR/data/variant_sets/$cohort.$context.collapsed_coding_csqs.tsv.gz \
       | gzip -c > $WRKDIR/data/variant_sets/$cohort.$context.recurrently_mutated_exons.tsv.gz"
  done
done


### Generate somatic and germline variant sets for burden testing
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
  for context in somatic germline; do
    case $context in
      germline)
        subset="RAS_loci"
        ;;
      somatic)
        subset="somatic_variants"
        ;;
    esac
    for suf in err log; do
      logfile=$WRKDIR/LSF/logs/generate_variant_sets_${cohort}_$context.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q normal -sla miket_sc \
      -J generate_variant_sets_${cohort}_$context \
      -o $WRKDIR/LSF/logs/generate_variant_sets_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/generate_variant_sets_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/populate_variant_sets.py \
         --vcf $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz \
         --sets-json $CODEDIR/refs/variant_set_criteria.$context.json \
       | gzip -c > $WRKDIR/data/variant_sets/$cohort.$context.burden_sets.tsv.gz"
  done
done


### Compute AC and AF matrixes by cancer type for all variants & variant sets in each cohort
# 1a. Coding variants by protein consequence
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

# 1b. All other single variants after excluding coding variants above
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

# 2. Recurrently mutated codons
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
      logfile=$WRKDIR/LSF/logs/get_recurrent_codon_freqs_${cohort}_$context.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
      -J get_recurrent_codon_freqs_${cohort}_$context \
      -o $WRKDIR/LSF/logs/get_recurrent_codon_freqs_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/get_recurrent_codon_freqs_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/calc_variant_set_freqs.py \
         --sets-tsv $WRKDIR/data/variant_sets/$cohort.$context.recurrently_mutated_codons.tsv.gz \
         --dosage-tsv $COHORTDIR/data/$cohort.$subset.dosage.tsv.gz \
         --sample-metadata $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
         --max-an $max_an \
         --outfile $WRKDIR/data/variant_set_freqs/$cohort.$context.recurrently_mutated_codons.freq.tsv.gz"
  done
done

# 3. Recurrently mutated exons
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
      logfile=$WRKDIR/LSF/logs/get_recurrent_exon_freqs_${cohort}_$context.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
      -J get_recurrent_exon_freqs_${cohort}_$context \
      -o $WRKDIR/LSF/logs/get_recurrent_exon_freqs_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/get_recurrent_exon_freqs_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/calc_variant_set_freqs.py \
         --sets-tsv $WRKDIR/data/variant_sets/$cohort.$context.recurrently_mutated_exons.tsv.gz \
         --dosage-tsv $COHORTDIR/data/$cohort.$subset.dosage.tsv.gz \
         --sample-metadata $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
         --max-an $max_an \
         --outfile $WRKDIR/data/variant_set_freqs/$cohort.$context.recurrently_mutated_exons.freq.tsv.gz"
  done
done

# 4. All variant sets
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
      logfile=$WRKDIR/LSF/logs/get_burden_set_freqs_${cohort}_$context.$suf
      if [ -e $logfile ]; then rm $logfile; fi
    done
    bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
      -J get_burden_set_freqs_${cohort}_$context \
      -o $WRKDIR/LSF/logs/get_burden_set_freqs_${cohort}_$context.log \
      -e $WRKDIR/LSF/logs/get_burden_set_freqs_${cohort}_$context.err \
      "$CODEDIR/scripts/data_processing/calc_variant_set_freqs.py \
         --sets-tsv $WRKDIR/data/variant_sets/$cohort.$context.burden_sets.tsv.gz \
         --dosage-tsv $COHORTDIR/data/$cohort.$subset.dosage.tsv.gz \
         --sample-metadata $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
         --max-an $max_an \
         --outfile $WRKDIR/data/variant_set_freqs/$cohort.$context.burden_sets.freq.tsv.gz"
  done
done

# # 5. Intra-gene somatic comutation pairs
# # TODO: finish implementing intra-gene restriction here
# for cohort in TCGA PROFILE; do
#   case $cohort in
#     TCGA)
#       COHORTDIR=$TCGADIR
#       ;;
#     PROFILE)
#       COHORTDIR=$PROFILEDIR
#       ;;
#   esac
#   # By definition, both mutations must each appear at ≥1% frequency for the pair
#   # to appear at ≥1% frequency
#   # We can use this definition to dramatically reduce the search space for computing
#   # pairwise comutation frequencies
#   # This requires pre-computed frequency info (generated above)
#   while read chrom start end gene; do
#     # Step 1. Build a list of all candidate mutations to consider for comutation
#     tabix \
#       -R <( echo -e "$chrom\t$start\t$end" ) \
#       $COHORTDIR/data/$cohort.somatic_variants.anno.clean.vcf.gz \
#     | cut -f3 > $TMPDIR/ras_vids.$gene.list
#     for subset in coding other; do
#       $CODEDIR/scripts/data_processing/filter_freq_table.py \
#         --freq-tsv $WRKDIR/data/variant_set_freqs/$cohort.somatic.${subset}_variants.freq.tsv.gz  \
#         --min-freq 0.01 \
#       | cut -f1 | sed '1d'
#     done | sort | uniq \
#     | fgrep -wf - \
#       <( zcat $WRKDIR/data/variant_sets/$cohort.somatic.collapsed_coding_csqs.tsv.gz \
#               $WRKDIR/data/variant_sets/$cohort.somatic.other_single_variants.tsv.gz ) \
#     | fgrep -wf $TMPDIR/ras_vids.$gene.list \
#     | awk -v OFS="\t" '{ print $1, $NF }' | cat <( echo -e "set_id\tvids" ) - \
#     > $TMPDIR/$cohort.all_comut_candidates.$gene.tsv
#     # Step 2. Compute comutation frequency for all candidates
#     for suf in err log; do
#       logfile=$WRKDIR/LSF/logs/get_somatic_comutation_freqs_${cohort}_$gene.$suf
#       if [ -e $logfile ]; then rm $logfile; fi
#     done
#     bsub -q big-multi -sla miket_sc -R "rusage[mem=24000]" \
#       -J get_somatic_comutation_freqs_${cohort}_$gene \
#       -o $WRKDIR/LSF/logs/get_somatic_comutation_freqs_${cohort}_$gene.log \
#       -e $WRKDIR/LSF/logs/get_somatic_comutation_freqs_${cohort}_$gene.err \
#       "$CODEDIR/scripts/data_processing/calc_comutation_freqs.py \
#          --sets-tsv $TMPDIR/$cohort.all_comut_candidates.$gene.tsv \
#          --dosage-tsv $COHORTDIR/data/$cohort.somatic_variants.dosage.tsv.gz \
#          --sample-metadata $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
#          --max-an 1 \
#          --outfile $WRKDIR/data/variant_set_freqs/$cohort.somatic.gene_comutations.$gene.freq.tsv.gz"
#   done < <( zcat $WRKDIR/../refs/RAS_genes.bed.gz )
# done
# # Collapse results across genes per cohort once complete
# for cohort in TCGA PROFILE; do
#   zcat $WRKDIR/data/variant_set_freqs/$cohort.somatic.gene_comutations.*.freq.tsv.gz \
#   | grep -ve '^set_id' | sort -Vk1,1 \
#   | cat <( zcat $WRKDIR/data/variant_set_freqs/$cohort.somatic.gene_comutations.NRAS.freq.tsv.gz | head -n1 ) - \
#   | bgzip -c > $WRKDIR/data/variant_set_freqs/$cohort.somatic.gene_comutations.freq.tsv.gz
# done

