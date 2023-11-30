#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Curate data from DFCI-PROFILE cohort for RAS modifiers study

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export BASEDIR=/data/gusev/PROFILE
export GTDIR=$BASEDIR/2020_2022_combined/IMPUTE_HQ
export CLINDIR=$BASEDIR/CLINICAL
export WRKDIR=/data/gusev/USERS/rlc47/PROFILE
export CODEDIR=$WRKDIR/../code/oncomod
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/sample_info LSF LSF/scripts LSF/logs ../refs refs; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Curate clinical information for patients of interest
# Get list of all IDs present in either imputed SNP and/or LOHGIC coding VCFs
tabix -H $GTDIR/PROFILE_COMB.22.HQ.vcf.gz \
| fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
> $WRKDIR/data/sample_info/PROFILE.imputed_snp.vcf.samples.list
tabix -H $WRKDIR/LOHGIC/data/PROFILE.LOHGIC.predicted_germline_coding_variants.vcf.gz \
| fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
> $WRKDIR/data/sample_info/PROFILE.oncopanel_lohgic.vcf.samples.list
cat \
  $WRKDIR/data/sample_info/PROFILE.imputed_snp.vcf.samples.list \
  $WRKDIR/data/sample_info/PROFILE.oncopanel_lohgic.vcf.samples.list \
| sort -V | uniq \
> $WRKDIR/data/sample_info/PROFILE.vcf.samples.list
# Curate EHR and get list of patients from cancer types of interest
$CODEDIR/scripts/data_processing/preprocess_dfci_profile_ehr.py \
  --id-map-tsv $CLINDIR/PROFILE_MRN_BL_PANEL.PBP.tab \
  --genomic-csv $CLINDIR/OncDRS/ALL_2021_11/GENOMIC_SPECIMEN.csv.gz \
  --dx-csv $CLINDIR/OncDRS/ALL_2021_11/CANCER_DIAGNOSIS_CAREG.csv.gz \
  --ancestry-csv $CLINDIR/PROFILE_2022_ANCESTRY.csv.gz \
  --survival-csv $CLINDIR/OncDRS/ALL_2021_11/PT_INFO_STATUS_REGISTRATION.csv.gz \
  --allfit-tsv $WRKDIR/LOHGIC/AllFIT/PROFILE.AllFIT_purity_estimates.tsv \
  --out-prefix $WRKDIR/data/sample_info/PROFILE. \
  --vcf-ids $WRKDIR/data/sample_info/PROFILE.vcf.samples.list \
  --priority-ids $WRKDIR/data/sample_info/PROFILE.oncopanel_lohgic.vcf.samples.list


### Subset VCFs to patients of interest and RAS loci
# # Recompress chromosomes 1-3 with bgzip & reindex
# for contig in $( seq 1 3 ); do
#   bsub -R 'rusage[mem=8000]' -q long -J compress_index_$contig \
#     "cd /data/gusev/PROFILE/2020_2022_combined/IMPUTE_HQ; \
#      . /PHShome/rlc47/.bashrc; \
#      zcat PROFILE_COMB.$contig.HQ.vcf.gz | bgzip -c > PROFILE_COMB.$contig.HQ.vcf.bgz; \
#      bcftools index PROFILE_COMB.$contig.HQ.vcf.bgz"
# done
# Extract imputed SNPs for samples & loci of interest
for contig in $( seq 1 22 ); do
  cat << EOF > $WRKDIR/LSF/scripts/extract_${contig}_imputed_snps.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools view \
  --min-ac 1 \
  --samples-file $WRKDIR/data/sample_info/PROFILE.ALL.samples.list \
  --regions-file $CODEDIR/refs/RAS_loci.plus_pathway.plus_GWAS.GRCh37.bed.gz \
  $GTDIR/PROFILE_COMB.$contig.HQ.vcf.gz \
| bcftools annotate -x FORMAT/DS \
| bcftools +fill-tags - \
  -O z -o $WRKDIR/data/PROFILE.imputed_snps.$contig.noGQs.vcf.gz \
  -- -t 'FORMAT/GQ:1=int(".")'
$CODEDIR/scripts/data_processing/gp2gq.py \
  $WRKDIR/data/PROFILE.imputed_snps.$contig.noGQs.vcf.gz \
  $WRKDIR/data/PROFILE.imputed_snps.$contig.vcf.gz
tabix -p vcf -f $WRKDIR/data/PROFILE.imputed_snps.$contig.vcf.gz
rm $WRKDIR/data/PROFILE.imputed_snps.$contig.noGQs.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/extract_${contig}_imputed_snps.sh
  rm $WRKDIR/LSF/logs/extract_${contig}_imputed_snps.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -n 2 -J PROFILE_extract_${contig}_imputed_snps \
    -o $WRKDIR/LSF/logs/extract_${contig}_imputed_snps.log \
    -e $WRKDIR/LSF/logs/extract_${contig}_imputed_snps.err \
    $WRKDIR/LSF/scripts/extract_${contig}_imputed_snps.sh
done
# Extract inferred coding variants for samples & loci of interest
for contig in $( seq 1 22 ); do
  bcftools view \
    -O z -o $WRKDIR/data/PROFILE.oncopanel_lohgic.$contig.vcf.gz \
    --min-ac 1 \
    --samples-file $WRKDIR/data/sample_info/PROFILE.ALL.samples.list \
    --force-samples \
    --regions-file <( zcat $CODEDIR/refs/RAS_loci.plus_pathway.plus_GWAS.GRCh37.bed.gz \
                      | awk -v contig=$contig '{ if ($1==contig) print }' ) \
    $WRKDIR/LOHGIC/data/PROFILE.LOHGIC.predicted_germline_coding_variants.vcf.gz
  tabix -p vcf -f $WRKDIR/data/PROFILE.oncopanel_lohgic.$contig.vcf.gz
done
# Merge imputed SNPs and inferred coding variants
for contig in $( seq 1 22 ); do
  cat << EOF > $WRKDIR/LSF/scripts/merge_imputed_lohgic_RAS_loci.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
$CODEDIR/scripts/data_processing/merge_profile_snps_lohgic.py \
  --lohgic-vcf $WRKDIR/data/PROFILE.oncopanel_lohgic.$contig.vcf.gz \
  --imputed-vcf $WRKDIR/data/PROFILE.imputed_snps.$contig.vcf.gz \
  --header $WRKDIR/refs/simple_hg19_header.vcf.gz \
  --outfile $WRKDIR/data/PROFILE.RAS_loci.$contig.vcf.gz
tabix -p vcf -f $WRKDIR/data/PROFILE.RAS_loci.$contig.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/merge_imputed_lohgic_RAS_loci.$contig.sh
  rm $WRKDIR/LSF/logs/merge_imputed_lohgic_RAS_loci.$contig.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -J PROFILE_merge_imputed_lohgic_RAS_loci_$contig \
    -o $WRKDIR/LSF/logs/merge_imputed_lohgic_RAS_loci.$contig.log \
    -e $WRKDIR/LSF/logs/merge_imputed_lohgic_RAS_loci.$contig.err \
    $WRKDIR/LSF/scripts/merge_imputed_lohgic_RAS_loci.$contig.sh
done
# # Merge combined VCFs for each chromosome into a single VCF and index the merged VCF
for contig in $( seq 1 22 ); do
  echo $WRKDIR/data/PROFILE.RAS_loci.$contig.vcf.gz
done > $WRKDIR/data/PROFILE.combined_vcf_shards.list
bcftools concat \
  --file-list $WRKDIR/data/PROFILE.combined_vcf_shards.list \
  -O z -o $WRKDIR/data/PROFILE.RAS_loci.vcf.gz
tabix -p vcf -f $WRKDIR/data/PROFILE.RAS_loci.vcf.gz


### Curate somatic data for patients of interest
# Download Gencode v19 GTF for extracting gene coordinates
wget -O /dev/stdout \
 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz \
| gunzip -c | grep -ve '^#' | sort -Vk1,1 -k4,4n -k5,5n | bgzip -c \
> $WRKDIR/../refs/gencode.v19.annotation.gtf.gz
tabix -f -p gff $WRKDIR/../refs/gencode.v19.annotation.gtf.gz
# Write list of donors missing somatic data
$CODEDIR/scripts/data_processing/find_dfci_profile_missing_somatic.py \
  --mutation-csv $CLINDIR/OncDRS/ALL_2021_11/GENOMIC_MUTATION_RESULTS.csv.gz \
  --cna-csv $CLINDIR/OncDRS/ALL_2021_11/GENOMIC_CNV_RESULTS.csv.gz \
  --samples-list $WRKDIR/data/sample_info/PROFILE.ALL.samples.list \
  --id-map-tsv $CLINDIR/PROFILE_MRN_BL_PANEL.PBP.tab \
  --out-prefix $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic
# Curate somatic data
$CODEDIR/scripts/data_processing/preprocess_dfci_profile_somatic.py \
  --mutation-csv $CLINDIR/OncDRS/ALL_2021_11/GENOMIC_MUTATION_RESULTS.csv.gz \
  --cna-csv $CLINDIR/OncDRS/ALL_2021_11/GENOMIC_CNV_RESULTS.csv.gz \
  --samples-list $WRKDIR/data/sample_info/PROFILE.ALL.samples.list \
  --no-mutation-data $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic.mut.list \
  --no-cna-data $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic.cna.list \
  --id-map-tsv $CLINDIR/PROFILE_MRN_BL_PANEL.PBP.tab \
  --genes-gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --priority-genes <( echo "KRAS" ) \
  --ref-fasta $WRKDIR/refs/GRCh37.fa \
  --lohgic $WRKDIR/LOHGIC/data/PROFILE.LOHGIC.annotated.final_predictions.tsv.gz \
  --header $WRKDIR/../refs/simple_hg19_header.somatic.vcf.gz \
  --outfile $WRKDIR/data/PROFILE.somatic_variants.vcf.gz \
  --out-tsv $WRKDIR/data/PROFILE.somatic_variants.tsv.gz
tabix -p vcf -f $WRKDIR/data/PROFILE.somatic_variants.vcf.gz


# # Curate selected PRS for PROFILE samples
# $CODEDIR/scripts/data_processing/curate_profile_prs.R \
#   $BASEDIR/23andme/DFCI_PRS_scores.txt \
#   $WRKDIR/data/sample_info/PROFILE.ALL.samples.list \
#   $CODEDIR/refs/PROFILE_selected_PRS.tsv \
#   $WRKDIR/data/PROFILE.PRS.tsv
# gzip -f $WRKDIR/data/PROFILE.PRS.tsv

