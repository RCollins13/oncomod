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
for SUBDIR in LSF LSF/logs LSF/scripts; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done
for SUBDIR in phase3 phase3/common_SNPs phase3/common_SNPs/LD_pruned; do
  if ! [ -e $KGDIR/$SUBDIR ]; then
    mkdir $KGDIR/$SUBDIR
  fi
done


### Prepare ref panel of LD-pruned common SNP markers for downstream analyses
# 1. Download VCFs with genotypes from 1kG phase 3 
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
# 2. Filter each to common biallelic SNPs
for contig in $( seq 1 22 ); do
  cat <<EOF > $WRKDIR/LSF/scripts/filter_auton_vcf.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $KGDIR/phase3
bcftools view \
  --include 'AF > 0.01 & AF < 1' \
  --types snps \
  -m 2 -M 2 \
  $KGDIR/phase3/ALL.chr$contig.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
  -Oz -o $KGDIR/phase3/common_SNPs/1000G.ph3.common_snps.chr$contig.vcf.gz
tabix -p vcf -f $KGDIR/phase3/common_SNPs/1000G.ph3.common_snps.chr$contig.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/filter_auton_vcf.$contig.sh
  rm $WRKDIR/LSF/logs/filter_auton_vcf.$contig.*
  bsub \
    -q normal -sla miket_sc -J filter_auton_vcf_$contig \
    -o $WRKDIR/LSF/logs/filter_auton_vcf.$contig.log \
    -e $WRKDIR/LSF/logs/filter_auton_vcf.$contig.err \
    $WRKDIR/LSF/scripts/filter_auton_vcf.$contig.sh
done
# 3. LD prune all common SNPs and retain sites
for contig in $( seq 1 22 ); do
  cat <<EOF > $WRKDIR/LSF/scripts/prune_auton_vcf.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $KGDIR/phase3/common_SNPs
module load plink/1.90b3
plink \
  --vcf $KGDIR/phase3/common_SNPs/1000G.ph3.common_snps.chr$contig.vcf.gz \
  --indep-pairwise 100kb 5 0.2 \
  --threads 4 \
  --memory 16000 \
  --out $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.chr$contig.pruned
bcftools view \
  --drop-genotypes \
  --include "ID=@$KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.chr$contig.pruned.prune.in" \
  --no-version \
  -Oz -o $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.chr$contig.pruned.vcf.gz \
  $KGDIR/phase3/common_SNPs/1000G.ph3.common_snps.chr$contig.vcf.gz
tabix -p vcf -f $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.chr$contig.pruned.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/prune_auton_vcf.$contig.sh
  rm $WRKDIR/LSF/logs/prune_auton_vcf.$contig.*
  bsub \
    -q big -sla miket_sc -R "rusage[mem=16000]" -n 4 \
    -sla miket_sc -J prune_auton_vcf_$contig \
    -o $WRKDIR/LSF/logs/prune_auton_vcf.$contig.log \
    -e $WRKDIR/LSF/logs/prune_auton_vcf.$contig.err \
    $WRKDIR/LSF/scripts/prune_auton_vcf.$contig.sh
done
# 4. Merge common sites across all contigs
for contig in $( seq 1 22 ); do
  echo $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.chr$contig.pruned.vcf.gz
done > $KGDIR/phase3/common_SNPs/LD_pruned/pruned_vcf_shards.list
bcftools concat \
  --file-list $KGDIR/phase3/common_SNPs/LD_pruned/pruned_vcf_shards.list \
  -O z -o $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.pruned.vcf.gz
tabix -p vcf -f $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.pruned.vcf.gz
# 5. Copy common SNPs reference to GCP (necessary for HMF data wrangling)
gsutil -m cp \
  $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.pruned.vcf.gz* \
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
  --regions-file $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.pruned.vcf.gz \
  /data/gusev/PROFILE/2020_2022_combined/IMPUTE_HQ/PROFILE_COMB.$contig.HQ.vcf.gz \
| bcftools annotate --no-version -x FORMAT/GP,FORMAT/DS \
  -Oz -o $PROFILEDIR/data/PROFILE.1000G_common_snps.$contig.vcf.gz
tabix -p vcf -f $PROFILEDIR/data/PROFILE.1000G_common_snps.$contig.vcf.gz
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
  echo $PROFILEDIR/data/PROFILE.1000G_common_snps.$contig.vcf.gz
done > $PROFILEDIR/data/1000G_common_snps.shards.list
bcftools concat \
  --file-list $PROFILEDIR/data/1000G_common_snps.shards.list \
| bcftools norm \
  --atomize \
  --check-ref x \
  --fasta-ref $TCGADIR/refs/GRCh37.fa \
  --threads 4 \
  -O z -o $PROFILEDIR/data/PROFILE.1000G_common_snps.vcf.gz
tabix -p vcf -f $PROFILEDIR/data/PROFILE.1000G_common_snps.vcf.gz


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
    <( bcftools view --regions $contig $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.pruned.vcf.gz ) \
  -Oz -o $TCGADIR/data/TCGA.1000G_common_snps.array_typed.$contig.vcf.gz \
  $TCGADIR/data/TCGA.array_typed.vcf.gz

# Imputed arrays
bcftools view \
  --no-version \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.array_imputed.samples.list \
  --regions-file \
    <( bcftools view --regions $contig $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.pruned.vcf.gz ) \
  -Oz -o $TCGADIR/data/TCGA.1000G_common_snps.array_imputed.$contig.vcf.gz \
  /data/gusev/TCGA/GENOTYPES/IMPUTED/$contig.vcf.gz

# Fake header for exome
bcftools view \
  --no-version \
  --samples-file $TCGADIR/data/sample_info/TCGA.ALL.exome.samples.list \
  --header-only \
  -Oz -o $TCGADIR/data/TCGA.1000G_common_snps.exome_dummy_empty.$contig.vcf.gz \
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
  --exome-vcf $TCGADIR/data/TCGA.1000G_common_snps.exome_dummy_empty.$contig.vcf.gz \
  --array-typed-vcf $TCGADIR/data/TCGA.1000G_common_snps.array_typed.$contig.vcf.gz \
  --array-imputed-vcf $TCGADIR/data/TCGA.1000G_common_snps.array_imputed.$contig.vcf.gz \
  --ref-fasta $TCGADIR/refs/GRCh37.fa \
  --header $TCGADIR/refs/simple_hg19_header.vcf.gz \
  --outfile $TCGADIR/data/TCGA.1000G_common_snps.$contig.vcf.gz \
  --verbose
tabix -p vcf -f $TCGADIR/data/TCGA.1000G_common_snps.$contig.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/merge_TCGA_imputed_typed_common_SNPs.$contig.sh
  rm $WRKDIR/LSF/logs/merge_TCGA_imputed_typed_common_SNPs.$contig.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -J merge_TCGA_imputed_typed_common_SNPs_$contig \
    -o $WRKDIR/LSF/logs/merge_TCGA_imputed_typed_common_SNPs.$contig.log \
    -e $WRKDIR/LSF/logs/merge_TCGA_imputed_typed_common_SNPs.$contig.err \
    $WRKDIR/LSF/scripts/merge_TCGA_imputed_typed_common_SNPs.$contig.sh
done
# 3. Combine SNP VCFs across all chromosomes
for contig in $( seq 1 22 ); do
  echo $TCGADIR/data/TCGA.1000G_common_snps.$contig.vcf.gz
done > $TCGADIR/data/1000G_common_snps.shards.list
bcftools concat \
  --file-list $TCGADIR/data/1000G_common_snps.shards.list \
| bcftools norm \
  --atomize \
  --check-ref x \
  --fasta-ref $TCGADIR/refs/GRCh37.fa \
  --threads 4 \
  -O z -o $TCGADIR/data/TCGA.1000G_common_snps.vcf.gz
tabix -p vcf -f $TCGADIR/data/TCGA.1000G_common_snps.vcf.gz


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
# 2. Merge VCFs across all samples
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
| bcftools norm \
  --atomize \
  --check-ref x \
  --fasta-ref $TCGADIR/refs/GRCh37.fa \
  --threads 4 \
  -Oz -o $HMFDIR/data/HMF.1000G_common_snps.vcf.gz
tabix -p vcf -f $HMFDIR/data/HMF.1000G_common_snps.vcf.gz
EOF
chmod a+x $WRKDIR/LSF/scripts/HMF_merge_common_snp_vcfs.sh
rm $WRKDIR/LSF/logs/HMF_merge_common_snp_vcfs.*
bsub \
  -q big -sla miket_sc \
  -n 4 -R 'rusage[mem=32000]' -J HMF_merge_common_snp_vcfs \
  -o $WRKDIR/LSF/logs/HMF_merge_common_snp_vcfs.log \
  -e $WRKDIR/LSF/logs/HMF_merge_common_snp_vcfs.err \
  $WRKDIR/LSF/scripts/HMF_merge_common_snp_vcfs.sh


### Merge SNPs across all cohorts
# Note: currently imposing no call rate filter (can be added downstream)
# May want to revisit this
bcftools merge \
  --no-version \
  --threads 4 \
  $PROFILEDIR/data/PROFILE.1000G_common_snps.vcf.gz \
  $TCGADIR/data/TCGA.1000G_common_snps.vcf.gz \
  $HMFDIR/data/HMF.1000G_common_snps.vcf.gz \
| bcftools +fill-tags -- -t AN,AC,AF,F_MISSING \
| bcftools view \
  --include 'INFO/AF > 0.01 & INFO/AF < 1' \
  -Oz -o $WRKDIR/data/all_cohorts.1000G_common_snps.vcf.gz


### Compute PC loadings with plink
module load plink/1.90b3
plink \
  --threads 4 \
  --memory 16000 \
  --vcf $WRKDIR/data/all_cohorts.1000G_common_snps.vcf.gz \
  --pca \
  --out $WRKDIR/data/all_cohorts.1000G_common_snps
