#!/usr/bin/env bash

################################
#    EGFR Modifiers Project    #
################################

# Copyright (c) 2023-Present Ryan L. Collins and the Gusev/Van Allen Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Annotate all VCFs prior to analysis

# Note: intended to be executed on the MGB ERISOne cluster


##################
### Basic prep ###
##################

### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export WRKDIR=/data/gusev/USERS/rlc47/EGFR_modifier_analysis
export CODEDIR=$WRKDIR/../EGFR_modifier_analysis/code/ras_modifiers
export VEP_CACHE=$WRKDIR/../refs/vep_cache
export VEP_PLUGINS=$WRKDIR/../code/vep_plugins
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in LSF LSF/logs LSF/scripts; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### NOTE: EXPECTS FULL RAS MODIFIER DATA CURATION TO BE ALREADY COMPLETE 
### SEE VERSION OF THIS FILE IN THE MAIN BRANCH


######################
### VEP annotation ###
######################

### Write script for generic VEP function call with plugins and options
# export PERL5LIB=$VEP_PLUGINS/loftee:$PERL5LIB
cat <<EOF > $WRKDIR/LSF/scripts/run_VEP.sh
#!/usr/bin/env bash
set -eu -o pipefail

# Extract CNVs larger than 100bp to be spiked in after VEP
# (VEP will exclude these given --max_sv_size 100)
cnv_records="$TMPDIR/\$( basename \$1 | sed 's/vcf.gz/cnvs.vcf_records/g' )"
bcftools view \$1 --no-header \
  --include 'INFO/SVTYPE ~ "DEL\|DUP\|AMP"' \
  -O v -o \$cnv_records

# Annotate with VEP
$WRKDIR/../code/ensembl-vep/vep \
  --input_file \$1 \
  --format vcf \
  --output_file STDOUT \
  --vcf \
  --force_overwrite \
  --species homo_sapiens \
  --assembly GRCh37 \
  --max_sv_size 100 \
  --offline \
  --no_stats \
  --cache \
  --dir_cache $VEP_CACHE/ \
  --cache_version 108 \
  --buffer_size 2500 \
  --dir_plugins $VEP_PLUGINS/ \
  --fasta $TCGADIR/refs/GRCh37.fa \
  --minimal \
  --sift b \
  --polyphen b \
  --nearest gene \
  --distance 10000 \
  --numbers \
  --hgvs \
  --no_escape \
  --symbol \
  --canonical \
  --domains \
  --plugin UTRAnnotator,file=$VEP_CACHE/UTRAnnotator/uORF_5UTR_GRCh37_PUBLIC.txt \
  --plugin dbNSFP,$VEP_CACHE/dbNSFP4.3a_grch37.gz,FATHMM_score,MPC_score \
  --plugin CADD,$VEP_CACHE/CADD/whole_genome_SNVs.tsv.gz,$VEP_CACHE/CADD/InDels.tsv.gz \
  --plugin LoF,loftee_path:$VEP_PLUGINS/loftee,human_ancestor_fa:$VEP_CACHE/loftee/human_ancestor.fa.gz,conservation_file:$VEP_CACHE/loftee/phylocsf_gerp.sql \
  --plugin SpliceAI,snv=$VEP_CACHE/SpliceAI/spliceai_scores.raw.snv.hg37.vcf.gz,indel=$VEP_CACHE/SpliceAI/spliceai_scores.raw.indel.hg37.vcf.gz \
  --custom $VEP_CACHE/All_GRCh37_RS.bw,GERP,bigwig \
  --custom $VEP_CACHE/GRCh37.100way.phastCons.bw,phastCons,bigwig \
  --custom $VEP_CACHE/hg19.100way.phyloP100way.bw,phyloP,bigwig \
  --custom $VEP_CACHE/gnomad/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_POPMAX \
  --custom $VEP_CACHE/gnomad/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz,gnomADe,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_POPMAX \
  --custom $VEP_CACHE/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG \
  --custom $VEP_CACHE/cosmic/cmc.vcf.gz,COSMIC,vcf,exact,0,COSMIC_GENE,COSMIC_AA_CHANGE,COSMIC_GENE_TIER,COSMIC_MUT_SIG,COSMIC_MUT_FREQ \
  --custom $VEP_CACHE/gencode_promoters.bed.gz,promoters,bed,overlap,0 \
  --custom $VEP_CACHE/ABC_enhancers.bed.gz,ABC_enhancers,bed,overlap,0 \
  --custom $VEP_CACHE/CCR.bed.gz,CCR,bed,overlap,0 \
  --custom $VEP_CACHE/samocha_regional_constraint.bed.gz,ExAC_regional_constraint,bed,overlap,0 \
  --custom $VEP_CACHE/GTEx_eQTLs.vcf.gz,GTEx,vcf,exact,0,GTEx_eGene,GTEx_eQTL_beta,GTEx_eQTL_tissue \
| cat - \$cnv_records \
| bcftools sort -O z -o \$2
rm \$cnv_records
tabix -p vcf -f \$2
EOF
chmod a+x $WRKDIR/LSF/scripts/run_VEP.sh


### Annotate VCFs
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
  esac
  for subset in EGFR_loci; do
    for suf in err log; do
      if [ -e $WRKDIR/LSF/logs/VEP_${cohort}_$subset.$suf ]; then
        rm $WRKDIR/LSF/logs/VEP_${cohort}_$subset.$suf
      fi
    done
    bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
      -J VEP_${cohort}_$subset \
      -o $WRKDIR/LSF/logs/VEP_${cohort}_$subset.log \
      -e $WRKDIR/LSF/logs/VEP_${cohort}_$subset.err \
      "$WRKDIR/LSF/scripts/run_VEP.sh \
         $COHORTDIR/data/$cohort.$subset.vcf.gz \
         $COHORTDIR/data/$cohort.$subset.anno.vcf.gz"
  done
done


#########################
### Clean VEP outputs ###
#########################

### Write script for cleaning generic VEP output (generated above)
cat <<EOF > $WRKDIR/LSF/scripts/clean_VEP.sh
#!/usr/bin/env bash
set -eu -o pipefail

$CODEDIR/scripts/data_processing/cleanup_vep.py \
  --gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --transcript-info $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
  \$1 stdout \
| grep -ve "^##bcftools" | grep -ve "^##CADD_" | grep -ve "^##UTRAnnotator_" \
| grep -ve "^##LoF_" | grep -ve "^##SpliceAI_" | grep -ve "^##VEP-command-line" \
| bgzip -c \
> \$2
tabix -f -p vcf \$2
EOF
chmod a+x $WRKDIR/LSF/scripts/clean_VEP.sh


### Clean up VCFs
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
  esac
  for subset in EGFR_loci; do
    for suf in err log; do
      if [ -e $WRKDIR/LSF/logs/VEP_cleanup_${cohort}_$subset.$suf ]; then
        rm $WRKDIR/LSF/logs/VEP_cleanup_${cohort}_$subset.$suf
      fi
    done
    bsub -q big -R "rusage[mem=16000]" -sla miket_sc \
      -J VEP_cleanup_${cohort}_$subset \
      -o $WRKDIR/LSF/logs/VEP_cleanup_${cohort}_$subset.log \
      -e $WRKDIR/LSF/logs/VEP_cleanup_${cohort}_$subset.err \
      "$WRKDIR/LSF/scripts/clean_VEP.sh \
         $COHORTDIR/data/$cohort.$subset.anno.vcf.gz \
         $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz"
  done
done


##############################
### Annotate in-sample AFs ###
##############################

### Write script to annotate in-cohort allele frequencies for all variants by population and cancer type
cat <<EOF > $WRKDIR/LSF/scripts/annotate_AFs.sh
#!/usr/bin/env bash
set -eu -o pipefail

$CODEDIR/scripts/data_processing/annotate_AFs.py \
  --sample-metadata \$2 \
  --sample-id-column \$3 \
  \$1 \
  \$4
tabix -f -p vcf \$4
EOF
chmod a+x $WRKDIR/LSF/scripts/annotate_AFs.sh


### Annotate AFs for all VCFs
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      sample_field=DONOR_ID
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      sample_field=PBP
      ;;
  esac
  for subset in EGFR_loci; do
    for suf in err log; do
      if [ -e $WRKDIR/LSF/logs/annotate_AFs_${cohort}_$subset.$suf ]; then
        rm $WRKDIR/LSF/logs/annotate_AFs_${cohort}_$subset.$suf
      fi
    done
    bsub -q long -sla miket_sc \
      -J annotate_AFs_${cohort}_$subset \
      -o $WRKDIR/LSF/logs/annotate_AFs_${cohort}_$subset.log \
      -e $WRKDIR/LSF/logs/annotate_AFs_${cohort}_$subset.err \
      "$WRKDIR/LSF/scripts/annotate_AFs.sh \
         $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz \
         $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
         $sample_field \
         $COHORTDIR/data/$cohort.$subset.anno.clean.wAF.vcf.gz"
  done
done


######################################
### Build simple genotype matrixes ###
######################################
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      sample_field=DONOR_ID
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      sample_field=PBP
      ;;
  esac
  for subset in EGFR_loci; do
    for suf in err log; do
      if [ -e $WRKDIR/LSF/logs/make_dosage_${cohort}_$subset.$suf ]; then
        rm $WRKDIR/LSF/logs/make_dosage_${cohort}_$subset.$suf
      fi
    done
    bsub -q normal -sla miket_sc \
      -J make_dosage_${cohort}_$subset \
      -o $WRKDIR/LSF/logs/make_dosage_${cohort}_$subset.log \
      -e $WRKDIR/LSF/logs/make_dosage_${cohort}_$subset.err \
      "$CODEDIR/scripts/data_processing/vcf2dosage.py \
         $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz - \
       | gzip -c > $COHORTDIR/data/$cohort.$subset.dosage.tsv.gz"
  done
done


################################################################
### Define sets of samples lacking DFCI/BWH Tier 1 mutations ###
################################################################
# Extract lists of samples
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      sample_field="donors"
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      sample_field="samples"
      ;;
  esac
  $CODEDIR/scripts/data_processing/define_control_samples.py \
    --vcf $COHORTDIR/data/$cohort.somatic_variants.anno.clean.vcf.gz \
    --criteria $CODEDIR/refs/EGFR_control_sample_criteria.json \
    --regions $CODEDIR/refs/EGFR_loci.GRCh37.bed.gz \
    --eligible-samples $COHORTDIR/data/sample_info/$cohort.LUAD.$sample_field.list \
    --exclude-samples $COHORTDIR/data/sample_info/$cohort.ALL.$sample_field.missing_somatic.list \
    --outfile $COHORTDIR/data/sample_info/$cohort.LUAD.eligible_EGFR_controls.list

done
# Summarize as table
for cancer in PDAC CRAD LUAD SKCM; do
  echo $cancer
  for cohort in TCGA PROFILE; do
    case $cohort in
      TCGA)
        COHORTDIR=$TCGADIR
        sample_field="donors"
        ;;
      PROFILE)
        COHORTDIR=$PROFILEDIR
        sample_field="samples"
        ;;
    esac
    elig_samps=$COHORTDIR/data/sample_info/$cohort.$cancer.$sample_field.list
    fgrep -wf $elig_samps \
      $COHORTDIR/data/sample_info/$cohort.ALL.eligible_controls.list \
    | wc -l | addcom
  done | paste -s -
done | paste - -
