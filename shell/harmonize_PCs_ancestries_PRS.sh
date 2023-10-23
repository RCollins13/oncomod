#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Unify germline SNPs from all cohorts for ancestry/relatedness inference & PRS

# Note: intended to be executed on the MGB ERISOne cluster


##################
### Basic prep ###
##################

### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export HMFDIR=/data/gusev/USERS/rlc47/HMF
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export KGDIR=/data/gusev/1000G/
export CODEDIR=$WRKDIR/../code/oncomod
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in LSF LSF/logs LSF/scripts data/1000G_SNPs; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done
for SUBDIR in phase3 phase3/common_SNPs phase3/common_SNPs/LD_pruned; do
  if ! [ -e $KGDIR/$SUBDIR ]; then
    mkdir $KGDIR/$SUBDIR
  fi
done


### Prepare ref panel data from 1kG phase 3
for contig in $( seq 1 22 ); do
  cat <<EOF > $WRKDIR/LSF/scripts/dl_auton_vcf.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $KGDIR/phase3
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$contig.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$contig.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
EOF
  chmod a+x $WRKDIR/LSF/scripts/dl_auton_vcf.$contig.sh
  rm $WRKDIR/LSF/logs/dl_auton_vcf.$contig.*
  bsub \
    -q filemove -sla miket_sc -J dl_auton_vcf_$contig \
    -o $WRKDIR/LSF/logs/dl_auton_vcf.$contig.log \
    -e $WRKDIR/LSF/logs/dl_auton_vcf.$contig.err \
    $WRKDIR/LSF/scripts/dl_auton_vcf.$contig.sh
done


### Identify SNP markers that are well-typed across all cohorts
# 1. Identify well-typed biallelic SNPs in TCGA array data
bcftools view \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.array_typed.samples.list \
  $TCGADIR/data/TCGA.array_typed.vcf.gz \
| bcftools +fill-tags - -- -t F_MISSING \
| bcftools view \
  --no-version -G --include 'F_MISSING < 0.02' -m 2 -M 2 --type snps \
  -Oz -o $TCGADIR/data/TCGA.well_typed_snps.sites.vcf.gz
tabix -p vcf -f $TCGADIR/data/TCGA.well_typed_snps.sites.vcf.gz
# 2. Intersect with well-imputed PROFILE SNPS
for contig in $( seq 1 22 ); do
  cat << EOF > $WRKDIR/LSF/scripts/PROFILE_get_well_typed_snps.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $PROFILEDIR
bcftools view \
  --no-version \
  --samples-file $PROFILEDIR/data/sample_info/PROFILE.ALL.samples.list \
  --regions-file $TCGADIR/data/TCGA.well_typed_snps.sites.vcf.gz \
  /data/gusev/PROFILE/2020_2022_combined/IMPUTE_HQ/PROFILE_COMB.$contig.HQ.vcf.gz \
| bcftools +fill-tags - -- -t F_MISSING \
| bcftools view \
  --no-version -G --include 'F_MISSING < 0.02' -m 2 -M 2 --type snps \
  -Oz -o $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_subset.$contig.vcf.gz
tabix -p vcf -f $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_subset.$contig.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/PROFILE_get_well_typed_snps.$contig.sh
  rm $WRKDIR/LSF/logs/PROFILE_get_well_typed_snps.$contig.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -n 2 -J PROFILE_get_well_typed_snps_$contig \
    -o $WRKDIR/LSF/logs/PROFILE_get_well_typed_snps.$contig.log \
    -e $WRKDIR/LSF/logs/PROFILE_get_well_typed_snps.$contig.err \
    $WRKDIR/LSF/scripts/PROFILE_get_well_typed_snps.$contig.sh
done
for contig in $( seq 1 22 ); do
  echo $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_subset.$contig.vcf.gz
done > $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_subset.vcfs.list
bcftools concat \
  --file-list $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_subset.vcfs.list \
  -Oz -o $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_subset.vcf.gz
# 3. LD prune these sites based on 1kG data
for contig in $( seq 1 22 ); do
  cat <<EOF > $WRKDIR/LSF/scripts/prune_wellGTed_snps.1000G.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
module load plink/1.90b3
bcftools view \
  --regions-file $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_subset.vcf.gz \
  --types snps \
  -m 2 -M 2 \
  $KGDIR/phase3/ALL.chr$contig.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
  -Oz -o $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.chr$contig.vcf.gz
tabix -p vcf -f $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.chr$contig.vcf.gz
plink \
  --vcf $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.chr$contig.vcf.gz \
  --indep-pairwise 100kb 5 0.2 \
  --threads 4 \
  --memory 6000 \
  --out $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.chr$contig.pruned
EOF
  chmod a+x $WRKDIR/LSF/scripts/prune_wellGTed_snps.1000G.$contig.sh
  rm $WRKDIR/LSF/logs/prune_wellGTed_snps.1000G.$contig.*
  bsub \
    -q normal -sla miket_sc -J prune_wellGTed_snps.1000G_$contig \
    -o $WRKDIR/LSF/logs/prune_wellGTed_snps.1000G.$contig.log \
    -e $WRKDIR/LSF/logs/prune_wellGTed_snps.1000G.$contig.err \
    $WRKDIR/LSF/scripts/prune_wellGTed_snps.1000G.$contig.sh
done
# 4. Compile final list of sites as VCF for subsetting
for contig in $( seq 1 22 ); do
  cat $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.chr$contig.pruned.prune.in
done | sort -V | uniq \
> $WRKDIR/data/well_typed_snps.TCGA_PROFILE_subset.pruned.ids.list
bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
  $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_subset.vcf.gz \
| bcftools view \
  --no-version \
  --include "ID=@$WRKDIR/data/well_typed_snps.TCGA_PROFILE_subset.pruned.ids.list" \
  -Oz -o $WRKDIR/data/well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz
tabix -p vcf -f $WRKDIR/data/well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz
# 4. Copy SNPs reference to GCP (necessary for HMF data wrangling)
gsutil -m cp \
  $WRKDIR/data/well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz* \
  gs://fc-d2262497-7f5e-49b4-89ff-77d322023f02/


### Filter PROFILE SNPs
# 1. Extract marker SNPs for samples of interest per chromosome 
for contig in $( seq 1 22 ); do
  cat << EOF > $WRKDIR/LSF/scripts/PROFILE_extract_common_snps.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $PROFILEDIR
bcftools view \
  --no-version \
  --samples-file $PROFILEDIR/data/sample_info/PROFILE.ALL.samples.list \
  --regions-file $WRKDIR/data/well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz \
  /data/gusev/PROFILE/2020_2022_combined/IMPUTE_HQ/PROFILE_COMB.$contig.HQ.vcf.gz \
| bcftools annotate \
  --no-version \
  -x FORMAT/GP,FORMAT/DS \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
  -Oz -o $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_PROFILE_subset.$contig.vcf.gz
tabix -p vcf -f $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_PROFILE_subset.$contig.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/PROFILE_extract_common_snps.$contig.sh
  rm $WRKDIR/LSF/logs/PROFILE_extract_common_snps.$contig.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -n 2 -J PROFILE_extract_common_snps_$contig \
    -o $WRKDIR/LSF/logs/PROFILE_extract_common_snps.$contig.log \
    -e $WRKDIR/LSF/logs/PROFILE_extract_common_snps.$contig.err \
    $WRKDIR/LSF/scripts/PROFILE_extract_common_snps.$contig.sh
done
# 2. Combine marker SNPs across all chromosomes
for contig in $( seq 1 22 ); do
  echo $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_PROFILE_subset.$contig.vcf.gz
done > $PROFILEDIR/data/well_typed_snps.TCGA_PROFILE_subset.shards.list
cat << EOF > $WRKDIR/LSF/scripts/PROFILE_merge_common_snp_vcfs.sh
bcftools concat \
  --file-list $PROFILEDIR/data/well_typed_snps.TCGA_PROFILE_subset.shards.list \
| bcftools norm \
  --atomize \
  --check-ref x \
  --fasta-ref $TCGADIR/refs/GRCh37.fa \
  --threads 4 \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
  -O z -o $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_PROFILE_subset.vcf.gz
tabix -p vcf -f $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_PROFILE_subset.vcf.gz
EOF
chmod a+x $WRKDIR/LSF/scripts/PROFILE_merge_common_snp_vcfs.sh
rm $WRKDIR/LSF/logs/PROFILE_merge_common_snp_vcfs.*
bsub \
  -q big -sla miket_sc \
  -n 4 -R 'rusage[mem=16000]' -J PROFILE_merge_common_snp_vcfs \
  -o $WRKDIR/LSF/logs/PROFILE_merge_common_snp_vcfs.log \
  -e $WRKDIR/LSF/logs/PROFILE_merge_common_snp_vcfs.err \
  $WRKDIR/LSF/scripts/PROFILE_merge_common_snp_vcfs.sh


## Filter TCGA SNPs
# 1. Subset typed arrays and imputed arrays to 1kG common SNPs
for contig in $( seq 1 22 ); do
  cat << EOF > $WRKDIR/LSF/scripts/TCGA_subset_common_SNPs_by_tech.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
set -eu -o pipefail
cd $TCGADIR

# Genotyped arrays
bcftools view \
  --no-version \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.array_typed.samples.list \
  --regions-file \
    <( bcftools view --regions $contig $WRKDIR/data/well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz ) \
  -Oz -o $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.array_typed.$contig.vcf.gz \
  $TCGADIR/data/TCGA.array_typed.vcf.gz

# Imputed arrays
bcftools view \
  --no-version \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.array_imputed.samples.list \
  --regions-file \
    <( bcftools view --regions $contig $WRKDIR/data/well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz ) \
  -Oz -o $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.array_imputed.$contig.vcf.gz \
  /data/gusev/TCGA/GENOTYPES/IMPUTED/$contig.vcf.gz

# Fake header for exome
bcftools view \
  --no-version \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.exome.samples.list \
  --header-only \
  -Oz -o $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.exome_dummy_empty.$contig.vcf.gz \
  $TCGADIR/data/TCGA.RAS_loci.exome.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/TCGA_subset_common_SNPs_by_tech.$contig.sh
  rm $WRKDIR/LSF/logs/TCGA_subset_common_SNPs_by_tech.$contig.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -sla miket_sc -J TCGA_subset_common_SNPs_by_tech_$contig \
    -o $WRKDIR/LSF/logs/TCGA_subset_common_SNPs_by_tech.$contig.log \
    -e $WRKDIR/LSF/logs/TCGA_subset_common_SNPs_by_tech.$contig.err \
    $WRKDIR/LSF/scripts/TCGA_subset_common_SNPs_by_tech.$contig.sh
done 
# 2. Unify common SNPs per chromosome
for contig in $( seq 1 22 ); do
  cat << EOF > $WRKDIR/LSF/scripts/merge_TCGA_imputed_typed_common_SNPs.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $TCGADIR
$CODEDIR/scripts/data_processing/merge_tcga_arrays_exomes.py \
  --sample-id-map $TCGADIR/data/sample_info/TCGA.ALL.id_map.tsv.gz \
  --exome-vcf $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.exome_dummy_empty.$contig.vcf.gz \
  --array-typed-vcf $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.array_typed.$contig.vcf.gz \
  --array-imputed-vcf $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.array_imputed.$contig.vcf.gz \
  --ref-fasta $TCGADIR/refs/GRCh37.fa \
  --header $TCGADIR/refs/simple_hg19_header.vcf.gz \
  --outfile $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.$contig.vcf.gz \
  --verbose
tabix -p vcf -f $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.$contig.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/merge_TCGA_imputed_typed_common_SNPs.$contig.sh
  rm $WRKDIR/LSF/logs/merge_TCGA_imputed_typed_common_SNPs.$contig.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -J merge_TCGA_imputed_typed_common_SNPs_$contig \
    -o $WRKDIR/LSF/logs/merge_TCGA_imputed_typed_common_SNPs.$contig.log \
    -e $WRKDIR/LSF/logs/merge_TCGA_imputed_typed_common_SNPs.$contig.err \
    $WRKDIR/LSF/scripts/merge_TCGA_imputed_typed_common_SNPs.$contig.sh
done
# 3. Combine SNP VCFs across all chromosomes and apply 99% call rate filter
for contig in $( seq 1 22 ); do
  echo $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.$contig.vcf.gz
done > $TCGADIR/data/well_typed_snps.TCGA_PROFILE_subset.shards.list
bcftools concat \
  --naive \
  --file-list $TCGADIR/data/well_typed_snps.TCGA_PROFILE_subset.shards.list \
| bcftools +fill-tags -- -t AN,AC,AF,F_MISSING \
| bcftools view \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.donors.list \
  -i 'INFO/F_MISSING < 0.02' \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
| bcftools norm \
  --atomize \
  --check-ref x \
  --fasta-ref $TCGADIR/refs/GRCh37.fa \
  --threads 4 \
  -O z -o $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.vcf.gz
tabix -p vcf -f $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.vcf.gz


### Filter HMF SNPs
# 1. Download per-sample VCFs subset in Terra
#    Note: This requires sample metadata manually downloaded from Terra
col_idxs=$( head -n1 $HMFDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
            | sed 's/\t/\n/g' | awk '{ if ($1 ~ "common_snps_vcf") print NR }' | paste -s -d, )
if ! [ -e $HMFDIR/data/sample_common_snp_vcfs ]; then
  mkdir $HMFDIR/data/sample_common_snp_vcfs
fi
sed '1d' $HMFDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
| cut -f$col_idxs | sed 's/\t/\n/g' \
| gsutil -m cp -I $HMFDIR/data/sample_common_snp_vcfs/
# 2. Merge VCFs across all samples & apply 98% call rate filter
find $HMFDIR/data/sample_common_snp_vcfs/ -name "*vcf.gz" \
> $HMFDIR/data/HMF.sample_common_snp_vcfs.list
cat << EOF > $WRKDIR/LSF/scripts/HMF_merge_common_snp_vcfs.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools merge \
  --file-list $HMFDIR/data/HMF.sample_common_snp_vcfs.list \
  --no-version \
  --missing-to-ref \
  --threads 4 \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
| bcftools view \
  --samples-file $HMFDIR/data/sample_info/HMF.ALL.samples.list \
| bcftools +fill-tags - -- -t F_MISSING \
| bcftools view \
  -i 'INFO/F_MISSING < 0.02' \
| bcftools norm \
  --atomize \
  --check-ref x \
  --fasta-ref $TCGADIR/refs/GRCh37.fa \
  --threads 4 \
  -Oz -o $HMFDIR/data/HMF.well_typed_snps.TCGA_PROFILE_subset.vcf.gz
tabix -p vcf -f $HMFDIR/data/HMF.well_typed_snps.TCGA_PROFILE_subset.vcf.gz
EOF
chmod a+x $WRKDIR/LSF/scripts/HMF_merge_common_snp_vcfs.sh
rm $WRKDIR/LSF/logs/HMF_merge_common_snp_vcfs.*
bsub \
  -q big -sla miket_sc \
  -n 4 -R 'rusage[mem=16000]' -J HMF_merge_common_snp_vcfs \
  -o $WRKDIR/LSF/logs/HMF_merge_common_snp_vcfs.log \
  -e $WRKDIR/LSF/logs/HMF_merge_common_snp_vcfs.err \
  $WRKDIR/LSF/scripts/HMF_merge_common_snp_vcfs.sh


### Merge 1kG genotypes for select loci for ground truth data in ancestry assignment
for contig in $( seq 1 22 ); do
  echo $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.chr$contig.vcf.gz
done > $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.shards.list
bcftools concat \
  --file-list $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.shards.list \
  -a --regions-file $WRKDIR/data/well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz \
| bcftools norm \
  --atomize \
  --check-ref x \
  --fasta-ref $TCGADIR/refs/GRCh37.fa \
  --threads 4 \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
  -O z -o $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz
tabix -p vcf -f $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz


### Merge SNPs across all cohorts
# 2. Merge genotypes for SNPs that are well-typed in all three cohorts
# Also add 1kG samples at this stage for ground truth ancestry inference data
bcftools merge \
  --no-version \
  --threads 4 \
  $PROFILEDIR/data/PROFILE.well_typed_snps.TCGA_PROFILE_subset.vcf.gz \
  $TCGADIR/data/TCGA.well_typed_snps.TCGA_PROFILE_subset.vcf.gz \
  $HMFDIR/data/HMF.well_typed_snps.TCGA_PROFILE_subset.vcf.gz \
  $WRKDIR/data/1000G_SNPs/1000G.ph3.well_typed_snps.TCGA_PROFILE_subset.pruned.vcf.gz \
| bcftools +fill-tags -- -t AN,AC,AF,F_MISSING \
| bcftools view \
  --include 'INFO/F_MISSING < 0.01' \
  -Oz -o $WRKDIR/data/all_cohorts.well_typed_snps.TCGA_PROFILE_subset.vcf.gz


### Compute PC loadings and kinship with plink2
cat << EOF > $WRKDIR/LSF/scripts/study_wide_PCA.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
module load plink/2.0a2.3
plink2 \
  --threads 8 \
  --memory 32000 \
  --vcf $WRKDIR/data/all_cohorts.well_typed_snps.TCGA_PROFILE_subset.vcf.gz \
  --geno 0.1 \
  --pca \
  --make-king-table \
  --king-table-filter 0.1 \
  --out $WRKDIR/data/all_cohorts.well_typed_snps.TCGA_PROFILE_subset
EOF
chmod a+x $WRKDIR/LSF/scripts/study_wide_PCA.sh
rm $WRKDIR/LSF/logs/study_wide_PCA.*
bsub \
  -q big-multi -sla miket_sc \
  -n 8 -R 'rusage[mem=32000]' -J study_wide_PCA \
  -o $WRKDIR/LSF/logs/study_wide_PCA.log \
  -e $WRKDIR/LSF/logs/study_wide_PCA.err \
  $WRKDIR/LSF/scripts/study_wide_PCA.sh

