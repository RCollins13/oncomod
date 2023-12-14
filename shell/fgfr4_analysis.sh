#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Focused investigation of germline FGFR4 association with KRAS tier 1 mutations

# Note: intended to be executed on the MGB ERISOne cluster


### Setup
# Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export HMFDIR=/data/gusev/USERS/rlc47/HMF
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export CODEDIR=$WRKDIR/../code/oncomod
cd $WRKDIR
# Prep directory structure
for dir in $WRKDIR/data/FGFR4 $WRKDIR/results/FGFR4 \
           $WRKDIR/plots/germline_somatic_assoc/FGFR4 \
           $WRKDIR/data/FGFR4/HMF_sample_vcfs; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done
# Ensure most recent version of OncoMod R package is installed from source
cd $CODEDIR && \
git pull && \
Rscript -e "install.packages('$CODEDIR/src/OncoModR_0.2.tar.gz', \
                             lib='~/R/x86_64-pc-linux-gnu-library/3.6', \
                             type='source', repos=NULL)" && \
cd -


### Compute single-variant association tests for all FGFR4 coding variants vs. all KRAS somatic endpoints
# Make input germline test sets
for cohort in TCGA PROFILE HMF; do
  zcat $WRKDIR/data/variant_sets/$cohort.germline.collapsed_coding_csqs.tsv.gz \
  | awk -v FS="\t" -v OFS="\t" '{ if ($2=="FGFR4") print $1, $NF }' \
  | sort -Vk1,1
done \
| $CODEDIR/scripts/germline_somatic_assoc/simple_collapse_variant_sets.py \
  --tsv-out $WRKDIR/data/FGFR4/FGFR4.coding_variants.sets.tsv
# Run pooled association tests for each FGFR4 coding variant
cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.sh
#!/usr/bin/env bash
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.pooled.R \
--sample-metadata $TCGADIR/data/sample_info/TCGA.ALL.sample_metadata.tsv.gz \
--somatic-ad $TCGADIR/data/TCGA.somatic_variants.dosage.tsv.gz \
--germline-ad $TCGADIR/data/TCGA.RAS_loci.dosage.tsv.gz \
--name TCGA \
--germline-gq $TCGADIR/data/TCGA.RAS_loci.GQ.tsv.gz \
--eligible-controls $TCGADIR/data/sample_info/TCGA.ALL.eligible_controls.list \
--sample-metadata $PROFILEDIR/data/sample_info/PROFILE.ALL.sample_metadata.tsv.gz \
--somatic-ad $PROFILEDIR/data/PROFILE.somatic_variants.dosage.tsv.gz \
--germline-ad $PROFILEDIR/data/PROFILE.RAS_loci.dosage.tsv.gz \
--name DFCI \
--germline-gq $PROFILEDIR/data/PROFILE.RAS_loci.GQ.tsv.gz \
--eligible-controls $PROFILEDIR/data/sample_info/PROFILE.ALL.eligible_controls.list \
--sample-metadata $HMFDIR/data/sample_info/HMF.ALL.sample_metadata.tsv.gz \
--somatic-ad $HMFDIR/data/HMF.somatic_variants.dosage.tsv.gz \
--germline-ad $HMFDIR/data/HMF.RAS_loci.dosage.tsv.gz \
--name HMF \
--germline-gq $HMFDIR/data/HMF.RAS_loci.GQ.tsv.gz \
--eligible-controls $HMFDIR/data/sample_info/HMF.ALL.eligible_controls.list \
--cancer-type CRAD \
--somatic-variant-sets $WRKDIR/data/variant_sets/test_sets/CRAD.KRAS.somatic_endpoints.tsv \
--germline-variant-sets $WRKDIR/data/FGFR4/FGFR4.coding_variants.sets.tsv \
--outfile $WRKDIR/results/FGFR4/pooled.CRAD.FGFR4_allelic_series.sumstats.tsv
gzip -f $WRKDIR/results/FGFR4/pooled.CRAD.FGFR4_allelic_series.sumstats.tsv
EOF
chmod a+x $WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.sh
for suf in err log; do
  logfile=$WRKDIR/LSF/logs/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.$suf
  if [ -e $logfile ]; then rm $logfile; fi
done
bsub -q big-multi -sla miket_sc -R "rusage[mem=24000]" -n 4 \
  -J germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series \
  -o $WRKDIR/LSF/logs/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.log \
  -e $WRKDIR/LSF/logs/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.err \
  "$WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.sh"


### Plot single-variant association results
$CODEDIR/scripts/plot/plot_FGFR4_results.R \
  $WRKDIR/results/FGFR4/pooled.CRAD.FGFR4_allelic_series.sumstats.tsv.gz \
  $WRKDIR/plots/germline_somatic_assoc/FGFR4


### Gather genotypes for all variants in the minimal FGFR4 haploblock 
### based on recombination data from LocusZoom
export fgfr4_haploblock="5:176507000-176593000"

## PROFILE
# Extract imputed SNPs in FGFR4 haploblock
bcftools view \
  --min-ac 1 \
  --samples-file $PROFILEDIR/data/sample_info/PROFILE.ALL.samples.list \
  --regions $fgfr4_haploblock \
  /data/gusev/PROFILE/2020_2022_combined/IMPUTE_HQ/PROFILE_COMB.5.HQ.vcf.gz \
| bcftools annotate -x FORMAT/DS \
| bcftools +fill-tags - \
  -O z -o $WRKDIR/data/FGFR4/PROFILE.imputed_snps.FGFR4.noGQs.vcf.gz \
  -- -t 'FORMAT/GQ:1=int(".")'
$CODEDIR/scripts/data_processing/gp2gq.py \
  $WRKDIR/data/FGFR4/PROFILE.imputed_snps.FGFR4.noGQs.vcf.gz \
  $WRKDIR/data/FGFR4/PROFILE.imputed_snps.FGFR4.wGQs.vcf.gz
bcftools annotate -x FORMAT/GP \
  -O z -o $WRKDIR/data/FGFR4/PROFILE.imputed_snps.FGFR4.vcf.gz \
  $WRKDIR/data/FGFR4/PROFILE.imputed_snps.FGFR4.wGQs.vcf.gz
tabix -p vcf -f $WRKDIR/data/FGFR4/PROFILE.imputed_snps.FGFR4.vcf.gz
rm \
  $WRKDIR/data/FGFR4/PROFILE.imputed_snps.FGFR4.noGQs.vcf.gz \
  $WRKDIR/data/FGFR4/PROFILE.imputed_snps.FGFR4.wGQs.vcf.gz
# Extract inferred FGFR4 coding variants
bcftools view \
  -O z -o $WRKDIR/data/FGFR4/PROFILE.oncopanel_lohgic.FGFR4.vcf.gz \
  --min-ac 1 \
  --samples-file $PROFILEDIR/data/sample_info/PROFILE.ALL.samples.list \
  --force-samples \
  --regions $fgfr4_haploblock \
  $PROFILEDIR/LOHGIC/data/PROFILE.LOHGIC.predicted_germline_coding_variants.vcf.gz
tabix -p vcf -f $WRKDIR/data/FGFR4/PROFILE.oncopanel_lohgic.FGFR4.vcf.gz
# Merge imputed SNPs and inferred coding variants
$CODEDIR/scripts/data_processing/merge_profile_snps_lohgic.py \
  --lohgic-vcf $WRKDIR/data/FGFR4/PROFILE.oncopanel_lohgic.FGFR4.vcf.gz \
  --imputed-vcf $WRKDIR/data/FGFR4/PROFILE.imputed_snps.FGFR4.vcf.gz \
  --header $PROFILEDIR/refs/simple_hg19_header.wGQ.vcf.gz \
  --outfile $WRKDIR/data/FGFR4/PROFILE.FGFR4.vcf.gz
tabix -p vcf -f $WRKDIR/data/FGFR4/PROFILE.FGFR4.vcf.gz

## TCGA
# Extract genotyped SNPs in FGFR4 haploblock
bcftools view \
  --min-ac 1  \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.array_typed.samples.list \
  --regions $fgfr4_haploblock \
  $TCGADIR/data/TCGA.array_typed.vcf.gz \
| bcftools +fill-tags - \
   -O z -o $WRKDIR/data/FGFR4/TCGA.array_typed.FGFR4.noGQs.vcf.gz \
   -- -t 'FORMAT/GQ:1=int(".")'
$CODEDIR/scripts/data_processing/fill_missing_vcf_format_values.py \
  --overwrite \
  $WRKDIR/data/FGFR4/TCGA.array_typed.FGFR4.noGQs.vcf.gz \
  $TCGADIR/data/sample_info/TCGA.ALL.typed_array_approx_GQ_per_sample.tsv.gz \
  $WRKDIR/data/FGFR4/TCGA.array_typed.FGFR4.vcf.gz
tabix -p vcf -f $WRKDIR/data/FGFR4/TCGA.array_typed.FGFR4.vcf.gz
rm $WRKDIR/data/FGFR4/TCGA.array_typed.FGFR4.noGQs.vcf.gz
# Extract imputed SNPs in FGFR4 haploblock
bcftools view \
  --min-ac 1 --exclude 'INFO/INFO < 0.8' \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.array_imputed.samples.list \
  --regions $fgfr4_haploblock \
  /data/gusev/TCGA/GENOTYPES/IMPUTED/5.vcf.gz \
| bcftools annotate -x FORMAT/ADS,FORMAT/DS \
| bcftools +fill-tags - \
  -Oz -o $WRKDIR/data/FGFR4/TCGA.array_imputed.FGFR4.noGQs.vcf.gz \
  -- -t 'FORMAT/GQ:1=int(".")'
$CODEDIR/scripts/data_processing/gp2gq.py \
  $WRKDIR/data/FGFR4/TCGA.array_imputed.FGFR4.noGQs.vcf.gz \
  $WRKDIR/data/FGFR4/TCGA.array_imputed.FGFR4.wGQs.vcf.gz
bcftools annotate -x FORMAT/GP \
  -O z -o $WRKDIR/data/FGFR4/TCGA.array_imputed.FGFR4.vcf.gz \
  $WRKDIR/data/FGFR4/TCGA.array_imputed.FGFR4.wGQs.vcf.gz
tabix -p vcf -f $WRKDIR/data/FGFR4/TCGA.array_imputed.FGFR4.vcf.gz
rm \
  $WRKDIR/data/FGFR4/TCGA.array_imputed.FGFR4.noGQs.vcf.gz \
  $WRKDIR/data/FGFR4/TCGA.array_imputed.FGFR4.wGQs.vcf.gz
# Use bcftools to stream WES data from gs:// bucket for FGFR4
gcloud auth application-default login
export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
cd $TCGADIR/data && \
bcftools view \
  --min-ac 1 \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.exome.samples.list \
  --regions $fgfr4_haploblock \
  gs://terra-workspace-archive-lifecycle/fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz \
| bcftools annotate -x FORMAT/AD,FORMAT/DP,FORMAT/PL,FORMAT/VAF \
  -O z -o $WRKDIR/data/FGFR4/TCGA.exome.FGFR4.noGQs.vcf.gz && \
cd -
# 3. Fill reference GQs per sample with sample-specific median homalt GQ
$CODEDIR/scripts/data_processing/fill_missing_vcf_format_values.py \
  --ignore-genotype \
  $WRKDIR/data/FGFR4/TCGA.exome.FGFR4.noGQs.vcf.gz \
  $TCGADIR/data/sample_info/TCGA.median_homalt_gqs.tsv.gz \
  $WRKDIR/data/FGFR4/TCGA.exome.FGFR4.vcf.gz
tabix -p vcf -f $WRKDIR/data/FGFR4/TCGA.exome.FGFR4.vcf.gz
rm $WRKDIR/data/FGFR4/TCGA.exome.FGFR4.noGQs.vcf.gz
# Merge VCFs for exomes and arrays for each chromosome
$CODEDIR/scripts/data_processing/merge_tcga_arrays_exomes.py \
  --sample-id-map $TCGADIR/data/sample_info/TCGA.ALL.id_map.tsv.gz \
  --exome-vcf $WRKDIR/data/FGFR4/TCGA.exome.FGFR4.vcf.gz \
  --array-typed-vcf $WRKDIR/data/FGFR4/TCGA.array_typed.FGFR4.vcf.gz \
  --array-imputed-vcf $WRKDIR/data/FGFR4/TCGA.array_imputed.FGFR4.vcf.gz \
  --ref-fasta $TCGADIR/refs/GRCh37.fa \
  --header $TCGADIR/refs/simple_hg19_header.wGQ.vcf.gz \
  --outfile $WRKDIR/data/FGFR4/TCGA.FGFR4.vcf.gz \
  --verbose
tabix -p vcf -f $WRKDIR/data/FGFR4/TCGA.FGFR4.vcf.gz

## HMF
# Download germline VCFs from Terra
# This requires sample metadata manually downloaded from Terra
col_idxs=$( head -n1 $HMFDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
            | sed 's/\t/\n/g' | awk '{ if ($1 ~ "fgfr4_vcf") print NR }' | paste -s -d, )
sed '1d' $HMFDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
| cut -f$col_idxs | sed 's/\t/\n/g' \
| gsutil -m cp -I $WRKDIR/data/FGFR4/HMF_sample_vcfs/
# Get list of all VCFs
find $WRKDIR/data/FGFR4/HMF_sample_vcfs/ -name "*vcf.gz" \
> $WRKDIR/data/FGFR4/HMF.sample_vcfs.list
# 1. Merge single-sample VCFs into cohort-wide VCFs per chromosome
bcftools merge \
  --file-list $WRKDIR/data/FGFR4/HMF.sample_vcfs.list \
  --no-version \
  --missing-to-ref \
  --threads 4 \
| bcftools view \
  --samples-file $HMFDIR/data/sample_info/HMF.ALL.samples.list \
| bcftools annotate \
  -x FORMAT/AD,FORMAT/DP,FORMAT/PL,FORMAT/MIN_DP,FORMAT/PGT,FORMAT/PID,FORMAT/RGQ,FORMAT/SB,INFO/AF,INFO/BaseQRankSum,INFO/ClippingRankSum,INFO/DB,INFO/DP,INFO/FS,INFO/MQ,INFO/MQRankSum,INFO/QD,INFO/ReadPosRankSum,INFO/SOR,INFO/ExcessHet \
| bcftools norm \
  --check-ref x \
  -m - \
  --fasta-ref $TCGADIR/refs/GRCh37.fa \
| bcftools +fill-tags - -- -t AN,AC,AF \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
| bcftools view \
  --trim-alt-alleles \
  -Oz -o $WRKDIR/data/FGFR4/HMF.FGFR4.missing_GQs.vcf.gz \
  --include 'AC > 0'
# 2. Fill reference GQs per sample with sample-specific median homalt GQ
$CODEDIR/scripts/data_processing/fill_missing_vcf_format_values.py \
  $WRKDIR/data/FGFR4/HMF.FGFR4.missing_GQs.vcf.gz \
  $HMFDIR/data/sample_info/HMF.median_homalt_gqs.tsv.gz \
  $WRKDIR/data/FGFR4/HMF.FGFR4.vcf.gz
tabix -p vcf -f $WRKDIR/data/FGFR4/HMF.FGFR4.vcf.gz
rm $WRKDIR/data/FGFR4/HMF.FGFR4.missing_GQs.vcf.gz*


## Merge FGFR4 genotypes across all cohorts
bcftools merge \
  -o $WRKDIR/data/FGFR4/all_cohorts.FGFR4.vcf.gz \
  -O z \
  $WRKDIR/data/FGFR4/TCGA.FGFR4.vcf.gz \
  $WRKDIR/data/FGFR4/PROFILE.FGFR4.vcf.gz \
  $WRKDIR/data/FGFR4/HMF.FGFR4.vcf.gz
tabix -p vcf -f $WRKDIR/data/FGFR4/all_cohorts.FGFR4.vcf.gz


## Impute genotypes with BEAGLE
# Download BEAGLE 5.4
mkdir $CODEDIR/../beagle
wget \
  --no-check-certificate \
  -P $CODEDIR/../beagle/ \
  https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
chmod a+x $CODEDIR/../beagle/beagle.22Jul22.46e.jar
# Download reference panel
wget \
  -P $CODEDIR/../beagle/ \
  https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.bref3/chr5.1kg.phase3.v5a.b37.bref3
# Impute genotypes for all samples
java -Xmx14g -jar $CODEDIR/../beagle/beagle.22Jul22.46e.jar \
  gt=$WRKDIR/data/FGFR4/all_cohorts.FGFR4.vcf.gz \
  map=$CODEDIR/../beagle/chr5.1kg.phase3.v5a.b37.bref3 \
  chrom=$fgfr4_haploblock \
  out=$WRKDIR/data/FGFR4/all_cohorts.FGFR4.imputed


## Make AD matrix based on 
$CODEDIR/scripts/data_processing/vcf2dosage.py \
  $WRKDIR/data/FGFR4/all_cohorts.FGFR4.vcf.gz - \
| gzip -c > $WRKDIR/data/FGFR4/all_cohorts.FGFR4.dosage.tsv.gz
