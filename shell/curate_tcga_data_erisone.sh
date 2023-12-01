#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Curate (some) data from TCGA cohort for RAS modifiers study

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export BASEDIR=/data/gusev/TCGA
export GTDIR=/data/gusev/TCGA/GENOTYPES
export WRKDIR=/data/gusev/USERS/rlc47/TCGA
export CODEDIR=$WRKDIR/../code/oncomod
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/sample_info LSF LSF/scripts LSF/logs refs misc; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Recode array-genotyped SNPs as VCF
module load plink/1.90b3
plink \
  --bfile $GTDIR/TCGA.NORMAL \
  --recode vcf bgz \
  --memory 16000 \
  --threads 4 \
  --out $WRKDIR/data/TCGA.array_typed
tabix -p vcf -f $WRKDIR/data/TCGA.array_typed.vcf.gz


### Define set of TCGA samples with complete data from cancer types of interest
# Get list of all IDs present in the WES VCF
gsutil -m cat \
  gs://fc-2cae157b-a485-483f-8f33-054640ba106f/2022_Final/101122_TCGA_10960_interval_process.tsv.gz \
| gunzip -c | cut -f2- | bgzip -c \
> $WRKDIR/data/101122_TCGA_10960_interval_process.tsv.gz
gsutil -m cp \
  gs://fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz.tbi \
  $WRKDIR/data/
gcloud auth application-default login
export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
cd $WRKDIR/data
tabix -H gs://fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz \
| fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
| fgrep -wf <( zcat $WRKDIR/data/101122_TCGA_10960_interval_process.tsv.gz \
               | head -n1 | cut -f2- | sed 's/\t/\n/g' ) \
> $WRKDIR/data/sample_info/TCGA.vcf.exome.samples.list
zcat $WRKDIR/refs/VanAllen.TCGA.WES_DeepVariant.sample_manifest.tsv.gz \
| fgrep -wf $WRKDIR/data/sample_info/TCGA.vcf.exome.samples.list \
| cut -f2 | cut -f1-3 -d\- | sort | uniq \
> $WRKDIR/data/sample_info/TCGA.vcf.exome.donors.list
# Get list of all IDs present in genotyped array VCF
tabix -H $WRKDIR/data//TCGA.array_typed.vcf.gz \
| fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
> $WRKDIR/data/sample_info/TCGA.vcf.array_typed.samples.gusev_IDs.list
# Get list of all IDs present in imputed VCF
tabix -H $GTDIR/IMPUTED/1.vcf.gz \
| fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
> $WRKDIR/data/sample_info/TCGA.vcf.array_imputed.samples.gusev_IDs.list
# Make map of Sasha's array IDs to canonical TCGA IDs
for tech in array_typed array_imputed; do
  cut --complement -f2 -d\- $WRKDIR/data/sample_info/TCGA.vcf.$tech.samples.gusev_IDs.list \
  > $WRKDIR/data/sample_info/TCGA.vcf.$tech.samples.list
  cut -f1-3 -d\- $WRKDIR/data/sample_info/TCGA.vcf.$tech.samples.list | sort | uniq \
  > $WRKDIR/data/sample_info/TCGA.vcf.$tech.donors.list
done
# Intersect samples with both exomes and genotyped arrays (not necessarily imputed)
$CODEDIR/scripts/data_processing/harmonize_tcga_samples.py \
  --exome-ids $WRKDIR/data/sample_info/TCGA.vcf.exome.samples.list \
  --exome-id-map $WRKDIR/refs/VanAllen.TCGA.WES_DeepVariant.sample_manifest.tsv.gz \
  --array-typed-ids $WRKDIR/data/sample_info/TCGA.vcf.array_typed.samples.gusev_IDs.list \
  --array-imputed-ids $WRKDIR/data/sample_info/TCGA.vcf.array_imputed.samples.gusev_IDs.list \
  --msi-tsv $WRKDIR/data/TCGA.MSI_mantis.Bonneville_2017.tsv.gz \
  --cdr-csv $BASEDIR/TCGA_CDR.csv \
  --tcga-tss-table $CODEDIR/refs/TCGA_TSS_codes.tsv.gz \
  --tcga-study-table $CODEDIR/refs/TCGA_study_codes.tsv.gz \
  --out-prefix $WRKDIR/data/sample_info/TCGA


### Curate clinical information for patients of interest
# Note: TCGA ancestry label assignments came from Carrot-Zhang et al., Cancer Cell, 2020
# https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020
# Filename: Broad_ancestry_PCA.txt
# This file has been renamed and relocated to $WRKDIR/data/TCGA.ancestry.tsv.gz
$CODEDIR/scripts/data_processing/preprocess_tcga_phenotypes.py \
  --id-map-tsv $WRKDIR/data/sample_info/TCGA.ALL.id_map.tsv.gz \
  --cdr-csv $BASEDIR/TCGA_CDR.csv \
  --tcga-study-table $CODEDIR/refs/TCGA_study_codes.tsv.gz \
  --ancestry-tsv $WRKDIR/data/TCGA.ancestry.tsv.gz \
  --pcs-txt $GTDIR/TCGA.COMBINED.QC.NORMAL.eigenvec \
  --purity-tsv $WRKDIR/data/TCGA.tumor_purity.Aran_2015.tsv.gz \
  --out-prefix $WRKDIR/data/sample_info/TCGA.


### Subset VCFs to patients of interest and RAS loci
# Estimate average GQ per sample from typed arrays based on the proportion of
# no-call SNPs per sample as a very rough proxy for sample genotyping quality
total_n_sites=$( bcftools query -f '%CHROM\n' $WRKDIR/data/TCGA.array_typed.vcf.gz | wc -l )
bcftools view \
  --samples-file $WRKDIR/data/sample_info/TCGA.ALL.$tech.samples.list \
  $WRKDIR/data/TCGA.array_typed.vcf.gz \
| bcftools query \
  -f '[%SAMPLE\n]' \
  --include 'GT=="mis"' \
| cat - $WRKDIR/data/sample_info/TCGA.ALL.$tech.samples.list \
| sort -V | uniq -c \
| awk -v OFS="\t" -v denom=$total_n_sites '{ print $2, $1/(denom+1) }' \
| awk -v OFS="\t" '{ gq=-10*log($2)/log(10); if(gq>99){gq=99}; print $1, gq }' \
| gzip -c \
> $WRKDIR/data/sample_info/TCGA.ALL.typed_array_approx_GQ_per_sample.tsv.gz
# Extract samples & loci of interest from genotyped arrays
export tech=array_typed
cat << EOF > $WRKDIR/LSF/scripts/extract_RAS_loci_variants_${tech}.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools view \
  --min-ac 1  \
  --samples-file $WRKDIR/data/sample_info/TCGA.ALL.$tech.samples.list \
  --regions-file $CODEDIR/refs/RAS_loci.plus_pathway.plus_GWAS.GRCh37.bed.gz \
  $WRKDIR/data/TCGA.array_typed.vcf.gz \
| bcftools +fill-tags - \
   -O z -o $WRKDIR/data/TCGA.RAS_loci.$tech.noGQs.vcf.gz \
   -- -t 'FORMAT/GQ:1=int(".")'
$CODEDIR/scripts/data_processing/fill_missing_vcf_format_values.py \
  --overwrite \
  $WRKDIR/data/TCGA.RAS_loci.$tech.noGQs.vcf.gz \
  $WRKDIR/data/sample_info/TCGA.ALL.typed_array_approx_GQ_per_sample.tsv.gz \
  $WRKDIR/data/TCGA.RAS_loci.$tech.vcf.gz
tabix -p vcf -f $WRKDIR/data/TCGA.RAS_loci.$tech.vcf.gz
EOF
chmod a+x $WRKDIR/LSF/scripts/extract_RAS_loci_variants_${tech}.sh
rm $WRKDIR/LSF/logs/extract_RAS_loci_variants_${tech}.*
bsub \
  -q short -R 'rusage[mem=6000]' -n 2 -J TCGA_extract_RAS_loci_variants_${tech} \
  -o $WRKDIR/LSF/logs/extract_RAS_loci_variants_${tech}.log \
  -e $WRKDIR/LSF/logs/extract_RAS_loci_variants_${tech}.err \
  $WRKDIR/LSF/scripts/extract_RAS_loci_variants_${tech}.sh
# Extract samples & loci of interest from imputed arrays
export tech=array_imputed
while read contig; do
  cat << EOF > $WRKDIR/LSF/scripts/extract_RAS_loci_variants_${tech}.${contig}.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools view \
  --min-ac 1 --exclude 'INFO/INFO < 0.8' \
  --samples-file $WRKDIR/data/sample_info/TCGA.ALL.$tech.samples.list \
  --regions-file $CODEDIR/refs/RAS_loci.plus_pathway.plus_GWAS.GRCh37.bed.gz \
  $GTDIR/IMPUTED/$contig.vcf.gz \
| bcftools annotate -x FORMAT/ADS,FORMAT/DS \
| bcftools +fill-tags - \
  -Oz -o $WRKDIR/data/TCGA.RAS_loci.$tech.$contig.noGQs.vcf.gz \
  -- -t 'FORMAT/GQ:1=int(".")'
$CODEDIR/scripts/data_processing/gp2gq.py \
  $WRKDIR/data/TCGA.RAS_loci.$tech.$contig.noGQs.vcf.gz \
  $WRKDIR/data/TCGA.RAS_loci.$tech.$contig.wGQs.vcf.gz
bcftools annotate -x FORMAT/GP \
  -O z -o $WRKDIR/data/TCGA.RAS_loci.$tech.$contig.vcf.gz \
  $WRKDIR/data/TCGA.RAS_loci.$tech.$contig.wGQs.vcf.gz
tabix -p vcf -f $WRKDIR/data/TCGA.RAS_loci.$tech.$contig.vcf.gz
rm \
  $WRKDIR/data/TCGA.RAS_loci.$tech.$contig.noGQs.vcf.gz \
  $WRKDIR/data/TCGA.RAS_loci.$tech.$contig.wGQs.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/extract_RAS_loci_variants_${tech}.${contig}.sh
  rm $WRKDIR/LSF/logs/extract_RAS_loci_variants_${tech}.${contig}.*
  bsub \
    -q short -R 'rusage[mem=6000]' -n 2 -J TCGA_extract_RAS_loci_variants_${tech}_${contig} \
    -o $WRKDIR/LSF/logs/extract_RAS_loci_variants_${tech}.${contig}.log \
    -e $WRKDIR/LSF/logs/extract_RAS_loci_variants_${tech}.${contig}.err \
    $WRKDIR/LSF/scripts/extract_RAS_loci_variants_${tech}.${contig}.sh
done < <( zcat $CODEDIR/refs/RAS_loci.plus_pathway.plus_GWAS.GRCh37.bed.gz \
          | fgrep -v "#" | cut -f1 | sort -V | uniq )
# Merge imputed array variants across all chromosomes
for contig in $( seq 1 22 ); do
  echo $WRKDIR/data/TCGA.RAS_loci.$tech.$contig.vcf.gz
done > $TMPDIR/$tech.contig_vcfs.list
bcftools concat \
  --file-list $TMPDIR/$tech.contig_vcfs.list \
  -O z -o $WRKDIR/data/TCGA.RAS_loci.$tech.vcf.gz
tabix -p vcf -f $WRKDIR/data/TCGA.RAS_loci.$tech.vcf.gz
# Define well-covered exome intervals within RAS loci of interest
$CODEDIR/scripts/data_processing/define_well_covered_targets.py \
  --coverage-matrix $WRKDIR/data/101122_TCGA_10960_interval_process.tsv.gz \
  --samples-list $WRKDIR/data/sample_info/TCGA.ALL.exome.samples.list \
  --min-frac-samples 0.9 \
  --min-frac-target 0.9 \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> $WRKDIR/refs/TCGA_WES.covered_intervals.bed.gz
tabix -p bed -f $WRKDIR/refs/TCGA_WES.covered_intervals.bed.gz
bedtools intersect -u \
  -a $WRKDIR/refs/TCGA_WES.covered_intervals.bed.gz \
  -b <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz \
             $CODEDIR/refs/NCI_RAS_pathway.genes.GRCh37.bed.gz \
        | fgrep -v "#" ) \
| bgzip -c > $WRKDIR/refs/TCGA_WES.covered_intervals.RAS_loci_plus_pathway.bed.gz
tabix -p bed -f $WRKDIR/refs/TCGA_WES.covered_intervals.RAS_loci_plus_pathway.bed.gz
# Use bcftools to stream WES data from gs:// bucket for samples & loci of interest
gsutil -m cp \
  gs://terra-workspace-archive-lifecycle/fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz.tbi \
  $WRKDIR/data/
gcloud auth application-default login
export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
cd $WRKDIR/data
bcftools view \
  --min-ac 1 \
  --samples-file $WRKDIR/data/sample_info/TCGA.ALL.exome.samples.list \
  --regions-file $WRKDIR/refs/TCGA_WES.covered_intervals.RAS_loci_plus_pathway.bed.gz \
  gs://terra-workspace-archive-lifecycle/fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz \
| bcftools annotate -x FORMAT/AD,FORMAT/DP,FORMAT/PL,FORMAT/VAF \
  -O z -o $WRKDIR/data/TCGA.RAS_loci.exome.noGQs.vcf.gz
# Get median GQ per sample for all autosomal biallelic homozygous SNPs
bcftools view \
  -m 2 -M 2 --types snps \
  $WRKDIR/data/TCGA.RAS_loci.exome.noGQs.vcf.gz \
| bcftools query -f '[%SAMPLE\t%GQ\n]' -i 'GT="AA"' \
| gzip -c \
> $WRKDIR/data/sample_info/TCGA.median_homalt_snp_gq.all_data.tsv.gz
$CODEDIR/utils/calc_sample_medians.py \
  $WRKDIR/data/sample_info/TCGA.median_homalt_snp_gq.all_data.tsv.gz \
| gzip -c > $WRKDIR/data/sample_info/TCGA.median_homalt_gqs.tsv.gz
# 3. Fill reference GQs per sample with sample-specific median homalt GQ
$CODEDIR/scripts/data_processing/fill_missing_vcf_format_values.py \
  --ignore-genotype \
  $WRKDIR/data/TCGA.RAS_loci.exome.noGQs.vcf.gz \
  $WRKDIR/data/sample_info/TCGA.median_homalt_gqs.tsv.gz \
  $WRKDIR/data/TCGA.RAS_loci.exome.vcf.gz
tabix -p vcf -f $WRKDIR/data/TCGA.RAS_loci.exome.vcf.gz
# Merge VCFs for exomes and arrays for each chromosome
for contig in $( seq 1 22 ); do
  cat << EOF > $WRKDIR/LSF/scripts/merge_arrays_exomes_RAS_loci.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
for tech in exome array_typed array_imputed; do
  bcftools view \
    --regions $contig \
    -Oz -o $WRKDIR/data/TCGA.RAS_loci.\$tech.$contig.vcf.gz \
    $WRKDIR/data/TCGA.RAS_loci.\$tech.vcf.gz
  tabix -p vcf -f $WRKDIR/data/TCGA.RAS_loci.\$tech.$contig.vcf.gz
done
$CODEDIR/scripts/data_processing/merge_tcga_arrays_exomes.py \
  --sample-id-map $WRKDIR/data/sample_info/TCGA.ALL.id_map.tsv.gz \
  --exome-vcf $WRKDIR/data/TCGA.RAS_loci.exome.$contig.vcf.gz \
  --array-typed-vcf $WRKDIR/data/TCGA.RAS_loci.array_typed.$contig.vcf.gz \
  --array-imputed-vcf $WRKDIR/data/TCGA.RAS_loci.array_imputed.$contig.vcf.gz \
  --ref-fasta $WRKDIR/refs/GRCh37.fa \
  --header $WRKDIR/refs/simple_hg19_header.wGQ.vcf.gz \
  --outfile $WRKDIR/data/TCGA.RAS_loci.$contig.vcf.gz \
  --verbose
tabix -p vcf -f $WRKDIR/data/TCGA.RAS_loci.$contig.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/merge_arrays_exomes_RAS_loci.$contig.sh
  rm $WRKDIR/LSF/logs/merge_arrays_exomes_RAS_loci.$contig.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -J TCGA_merge_arrays_exomes_RAS_loci.$contig \
    -o $WRKDIR/LSF/logs/merge_arrays_exomes_RAS_loci.$contig.log \
    -e $WRKDIR/LSF/logs/merge_arrays_exomes_RAS_loci.$contig.err \
    $WRKDIR/LSF/scripts/merge_arrays_exomes_RAS_loci.$contig.sh
done
# Merge tech-integrated germline variants across all chromosomes
for contig in $( seq 1 22 ); do
  echo $WRKDIR/data/TCGA.RAS_loci.$contig.vcf.gz
done > $WRKDIR/data/TCGA.RAS_loci.sharded_per_contig.vcfs.list
bcftools concat \
  --naive \
  --file-list $WRKDIR/data/TCGA.RAS_loci.sharded_per_contig.vcfs.list \
  -Oz -o $WRKDIR/data/TCGA.RAS_loci.vcf.gz
tabix -p vcf -f $WRKDIR/data/TCGA.RAS_loci.vcf.gz


### Curate somatic data for patients of interest
# Preprocess mc3 data by subsetting on samples of interest to improve pandas IO speed
cat \
  <( zcat $WRKDIR/data/mc3.v0.2.8.PUBLIC.maf.gz | head -n1 ) \
  <( zcat $WRKDIR/data/mc3.v0.2.8.PUBLIC.maf.gz \
     | fgrep -wf $WRKDIR/data/sample_info/TCGA.ALL.donors.list ) \
| gzip -c > $WRKDIR/data/mc3.v0.2.8.PUBLIC.PDAC_CRAD_LUAD_SKCM.maf.gz
# Pre-filter Affy SNP array CNA calls
# Note: must have been downloaded by cancer type from GDAC Firehose
# and uploaded to $WRKDIR/data/TCGA_CNA
cat $( find $WRKDIR/data/TCGA_CNA -name "*seg.seg.txt" | paste -s ) \
| fgrep -v "Num_Probes" | awk -v OFS="\t" \
  '{ cnv="NONE"; if ($6 <= -0.4150375) cnv="DEL"; if ($6 >= 0.3219281) cnv="AMP"; \
     if ($5>=10 && cnv!="NONE") print $2, $3, $4, cnv, $1 }' \
| sort -Vk1,1 -k2,2n -k3,3n -k5,5V -k4,4V \
| cat <( echo -e "#chrom\tstart\tend\tCNA\tsample" ) - | bgzip -c \
> $WRKDIR/data/TCGA.CNA.b19.bed.gz
# Write lists of donors with no somatic data available callset
zcat $WRKDIR/data/mc3.v0.2.8.PUBLIC.PDAC_CRAD_LUAD_SKCM.maf.gz \
| cut -f16 | cut -f1-3 -d\- | sort | uniq \
| fgrep -wvf - $WRKDIR/data/sample_info/TCGA.ALL.donors.list \
> $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.mc3.list
cat $( find $WRKDIR/data/TCGA_CNA -name "*seg.seg.txt" | paste -s ) \
| fgrep -v "Num_Probes" | cut -f1-3 -d\- | sort | uniq \
| fgrep -wvf - $WRKDIR/data/sample_info/TCGA.ALL.donors.list \
> $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.cna.list
cat \
  $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.mc3.list \
  $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.cna.list \
| sort | uniq \
> $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.list
# Curate somatic data
$CODEDIR/scripts/data_processing/preprocess_tcga_somatic.py \
  --mc3-tsv $WRKDIR/data/mc3.v0.2.8.PUBLIC.PDAC_CRAD_LUAD_SKCM.maf.gz \
  --cna-bed $WRKDIR/data/TCGA.CNA.b19.bed.gz \
  --donors-list $WRKDIR/data/sample_info/TCGA.ALL.donors.list \
  --no-mutation-data $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.mc3.list \
  --no-cna-data $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.cna.list \
  --genes-gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --priority-genes <( echo "KRAS" ) \
  --ref-fasta $WRKDIR/refs/GRCh37.fa \
  --header $WRKDIR/../refs/simple_hg19_header.somatic.vcf.gz \
  --outfile $WRKDIR/data/TCGA.somatic_variants.vcf.gz \
  --out-tsv $WRKDIR/data/TCGA.somatic_variants.tsv.gz
tabix -p vcf -f $WRKDIR/data/TCGA.somatic_variants.vcf.gz


# Summarize somatic variant status by gene & cancer type
for cancer in PDAC CRAD LUAD; do
  n_samp=$( fgrep -wvf \
              $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.list \
              ${WRKDIR}/data/sample_info/TCGA.$cancer.donors.list | wc -l )
  while read gene; do
    zcat $WRKDIR/data/TCGA.somatic_variants.tsv.gz \
    | fgrep -wf ${WRKDIR}/data/sample_info/TCGA.$cancer.donors.list \
    | fgrep -wvf $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.list \
    | awk -v gene=$gene -v FS="\t" '{ if ($9==gene && ($15~/AMP|DEL|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/)) print $7 }' \
    | sort | uniq | wc -l | awk -v n=$n_samp '{ print $1/n }'
  done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 )
  zcat $WRKDIR/data/TCGA.somatic_variants.tsv.gz \
    | fgrep -wf ${WRKDIR}/data/sample_info/TCGA.$cancer.donors.list \
    | fgrep -wvf $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.list \
    | awk -v FS="\t" '{ if ($9~/NRAS|HRAS|KRAS/ && $15~/AMP|DEL|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/) print $7 }' \
    | sort | uniq | wc -l | awk -v n=$n_samp '{ print $1/n }'
done | paste - - - -

