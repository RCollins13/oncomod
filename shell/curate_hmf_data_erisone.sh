#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Curate data from HMF cohort for RAS modifiers study

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export BASEDIR=/data/gusev/PROFILE
export WRKDIR=/data/gusev/USERS/rlc47/HMF
export CODEDIR=$WRKDIR/../code/oncomod
export KGDIR=/data/gusev/1000G/
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/sample_vcfs data/sample_somatic_vcf_subsets \
              data/sample_info data/all_somatic data/ancestry_inference \
              LSF LSF/scripts LSF/logs ../refs refs; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Collect Purple somatic QC (including sex inference & MSI status)
# Note that this has to be performed prior to any downstream data curation
# 1. Download & extract HMF somatic data
gsutil -m cp \
  gs://hmf-dr-355-us-central1/somatic.tar.gz \
  $WRKDIR/data/all_somatic/
tar -xzvf $WRKDIR/data/all_somatic/somatic.tar.gz
rm $WRKDIR/data/all_somatic/somatic.tar.gz
# 2. Interate over all samples and combine PURPLE somatic QC
find $WRKDIR/data/all_somatic/ -name "*.purple.purity.tsv" \
| head -n1 | xargs -I {} cat {} | head -n1 | sed 's/#//g' \
| awk -v OFS="\t" '{ print "#sample", $0 }' \
> $WRKDIR/data/all_somatic/HMF.all.purple.somatic_qc.tsv
while read file; do
  sid=$( basename $file | sed 's/.purple.purity.tsv//g' )
  cat $file | tail -n1 | paste <( echo -e "$sid" ) -
done < <( find $WRKDIR/data/all_somatic/ -name "*.purple.purity.tsv" ) \
| sort -Vk1,1 >> $WRKDIR/data/all_somatic/HMF.all.purple.somatic_qc.tsv
gzip -f $WRKDIR/data/all_somatic/HMF.all.purple.somatic_qc.tsv


### Download germline VCFs from Terra
# This requires sample metadata manually downloaded from Terra
col_idxs=$( head -n1 $WRKDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
            | sed 's/\t/\n/g' | awk '{ if ($1 ~ "ras_loci_vcf") print NR }' | paste -s -d, )
sed '1d' $WRKDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
| cut -f$col_idxs | sed 's/\t/\n/g' \
| gsutil -m cp -I $WRKDIR/data/sample_vcfs/
# Get list of all VCFs
find $WRKDIR/data/sample_vcfs/ -name "*vcf.gz" > $WRKDIR/data/HMF.sample_vcfs.list
# Get list of samples from all VCFs
while read vcf; do
  bcftools query -l $vcf
done < $WRKDIR/data/HMF.sample_vcfs.list | sort -V | uniq \
> $WRKDIR/data/HMF.ss_germline_vcfs.all_samples.list


### Curate clinical information for patients of interest
# 1. Download HMF metadata
gsutil -m cp \
  gs://hmf-dr-355-us-central1/metadata.tsv \
  gs://hmf-dr-355-us-central1/manifest.json \
  $WRKDIR/data/sample_info/
# 2. Curate metadata for patients with relevant cancer types and MSS tumors
$CODEDIR/scripts/data_processing/preprocess_hmf_phenotypes.py \
  --metadata $WRKDIR/data/sample_info/metadata.tsv \
  --purple-qc $WRKDIR/data/all_somatic/HMF.all.purple.somatic_qc.tsv.gz \
  --vcf-ids $WRKDIR/data/HMF.ss_germline_vcfs.all_samples.list \
  --hmf-json $WRKDIR/data/sample_info/manifest.json \
  --out-prefix $WRKDIR/data/sample_info/HMF.


### Curate somatic variants in regions of interest
# Note: assumes somatic data has been downloaded (as above)
# 1. Subset each sample's Purple output to KRAS
while read sid; do
  tid=$( echo $sid | sed 's/R$/T/g' \
         | fgrep -f - $WRKDIR/data/sample_info/metadata.tsv \
         | head -n1 | cut -f1)
  if [ -z $tid ]; then
    echo -e "$sid matched tumor cannot be found?"
  else
    cat << EOF > $WRKDIR/data/all_somatic/$tid/get_KRAS_somatic_variants.$sid.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools view \
  --samples $tid \
  --include 'FILTER == "PASS"' \
  --regions $( zcat $WRKDIR/../refs/RAS_genes.bed.gz | fgrep -w KRAS \
               | awk '{ print $1":"$2-10000"-"$3+10000 }' ) \
  $WRKDIR/data/all_somatic/$tid/purple/$tid.purple.somatic.vcf.gz \
| bcftools norm \
  -m - \
  --check-ref x \
  --fasta-ref /data/gusev/USERS/rlc47/TCGA/refs/GRCh37.fa \
| bcftools reheader \
  -s <( echo "$sid" ) \
| bcftools annotate \
  -x FORMAT/AD,FORMAT/AF,FORMAT/DP,FORMAT/RABQ,FORMAT/RAD,FORMAT/RC_CNT,FORMAT/RC_IPC,FORMAT/RC_JIT,FORMAT/RC_QUAL,FORMAT/RDP,FORMAT/SB,INFO/BIALLELIC,INFO/GND_FREQ,INFO/HOTSPOT,INFO/IMPACT,INFO/KT,INFO/LPS,INFO/LPS_RC,INFO/MAPPABILITY,INFO/MH,INFO/MSG,INFO/NEAR_HOTSPOT,INFO/PAVE_TI,INFO/PON_COUNT,INFO/PON_MAX,INFO/PURPLE_AF,INFO/PURPLE_CN,INFO/PURPLE_GERMLINE,INFO/PURPLE_MACN,INFO/PURPLE_VCN,INFO/RC,INFO/RC_IDX,INFO/RC_LF,INFO/RC_MH,INFO/RC_NM,INFO/RC_REPC,INFO/RC_REPS,INFO/RC_RF,INFO/REPORTED,INFO/REP_C,INFO/REP_S,INFO/SUBCL,INFO/TIER,INFO/TNC \
  -Oz -o $WRKDIR/data/sample_somatic_vcf_subsets/$sid.KRAS.simple.somatic.vcf.gz
tabix -p vcf -f $WRKDIR/data/sample_somatic_vcf_subsets/$sid.KRAS.simple.somatic.vcf.gz
EOF
    chmod a+x $WRKDIR/data/all_somatic/$tid/get_KRAS_somatic_variants.$sid.sh
    rm \
      $WRKDIR/data/all_somatic/$tid/get_KRAS_somatic_variants.$sid.log \
      $WRKDIR/data/all_somatic/$tid/get_KRAS_somatic_variants.$sid.err
    bsub \
      -q vshort -sla miket_sc -J get_KRAS_variants_$sid \
      -o $WRKDIR/data/all_somatic/$tid/get_KRAS_somatic_variants.$sid.log \
      -e $WRKDIR/data/all_somatic/$tid/get_KRAS_somatic_variants.$sid.err \
      $WRKDIR/data/all_somatic/$tid/get_KRAS_somatic_variants.$sid.sh
  fi
done < $WRKDIR/data/sample_info/HMF.ALL.samples.list
# 2. Merge somatic KRAS VCFs across all samples
cat << EOF > $WRKDIR/LSF/scripts/merge_somatic_KRAS_variants.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
find $WRKDIR/data/sample_somatic_vcf_subsets/ -name "*vcf.gz" \
> $WRKDIR/data/sample_somatic_vcf_subsets/all_somatic_KRAS_vcfs.list
bcftools merge \
  --file-list $WRKDIR/data/sample_somatic_vcf_subsets/all_somatic_KRAS_vcfs.list \
  --no-version \
  --missing-to-ref \
  --threads 4 \
| bcftools norm \
  -m - \
  --check-ref x \
  --fasta-ref /data/gusev/USERS/rlc47/TCGA/refs/GRCh37.fa \
| bcftools view \
  --samples-file $WRKDIR/data/sample_info/HMF.ALL.samples.list \
| bcftools +fill-tags - -- -t AN,AC,AF \
| bcftools annotate \
  -x FORMAT/SC_INS_COUNT,FORMAT/ABQ,FORMAT/UMI_CNT \
| bcftools view \
  -Oz -o $WRKDIR/data/HMF.somatic_variants.snvs_indels.vcf.gz \
  --include 'AC > 0'
tabix -p vcf -f $WRKDIR/data/HMF.somatic_variants.snvs_indels.vcf.gz
EOF
chmod a+x $WRKDIR/LSF/scripts/merge_somatic_KRAS_variants.sh
rm $WRKDIR/LSF/logs/merge_somatic_KRAS_variants.*
bsub \
  -q big -sla miket_sc \
  -n 4 -R 'rusage[mem=16000]' -J merge_somatic_KRAS_variants \
  -o $WRKDIR/LSF/logs/merge_somatic_KRAS_variants.log \
  -e $WRKDIR/LSF/logs/merge_somatic_KRAS_variants.err \
  $WRKDIR/LSF/scripts/merge_somatic_KRAS_variants.sh
# 3. Collect list of samples with KRAS CNAs
# Note: for our purposes, we consider a tumor as KRAS-amplified/-deleted if
# the KRAS copy number is at least 50% greater or less than the overall tumor ploidy.
# This is not perfect but likely good enough as a first pass.
while read sid; do
  tid=$( echo $sid | sed 's/R$/T/g' \
         | fgrep -f - $WRKDIR/data/sample_info/metadata.tsv \
         | head -n1 | cut -f1 )
  if [ -z $tid ]; then
    echo -e "$sid matched tumor cannot be found?"
    continue
  fi
  ploidy=$( cut -f5 $WRKDIR/data/all_somatic/$tid/purple/$tid.purple.purity.tsv | tail -n1 )
  awk -v FS="\t" '{ if ($4=="KRAS") print }' \
    $WRKDIR/data/all_somatic/$tid/purple/$tid.purple.cnv.gene.tsv \
  | awk -v FS="\t" -v OFS="\t" -v ploidy=$ploidy -v sid=$sid \
    '{ if ($6>=(1.5*ploidy)) print sid, "KRAS", "AMP"; else if ($5<=(0.5*ploidy)) print sid, "KRAS", "DEL" }'
done < $WRKDIR/data/sample_info/HMF.ALL.samples.list \
> $WRKDIR/data/HMF.KRAS_CNA_carriers.tsv
# 4. Unify somatic SNVs/indels and CNAs in KRAS
$CODEDIR/scripts/data_processing/preprocess_hmf_somatic.py \
  --muts-vcf $WRKDIR/data/HMF.somatic_variants.snvs_indels.vcf.gz \
  --cna-tsv $WRKDIR/data/HMF.KRAS_CNA_carriers.tsv \
  --genes-gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --priority-genes <( echo "KRAS" ) \
  --ref-fasta /data/gusev/USERS/rlc47/TCGA/refs/GRCh37.fa \
  --header $WRKDIR/../refs/simple_hg19_header.somatic.vcf.gz \
  --outfile $WRKDIR/data/HMF.somatic_variants.vcf.gz
tabix -p vcf -f $WRKDIR/data/HMF.somatic_variants.vcf.gz
# Since all HMF samples have somatic information, we need to make an empty list
# of samples missing somatic info (as this is expected by some downstream steps)
touch $WRKDIR/data/sample_info/HMF.ALL.samples.missing_somatic.list


### Curate germline variants in regions of interest
# 1. Merge single-sample VCFs into cohort-wide VCFs per chromosome
for contig in $( seq 1 22 ) X Y; do
  cat << EOF > $WRKDIR/LSF/scripts/merge_sample_vcfs.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools merge \
  --file-list $WRKDIR/data/HMF.sample_vcfs.list \
  --regions $contig \
  --no-version \
  --missing-to-ref \
  --threads 4 \
| bcftools view \
  --samples-file $WRKDIR/data/sample_info/HMF.ALL.samples.list \
| bcftools annotate \
  -x FORMAT/AD,FORMAT/DP,FORMAT/PL,FORMAT/MIN_DP,FORMAT/PGT,FORMAT/PID,FORMAT/RGQ,FORMAT/SB,INFO/AF,INFO/BaseQRankSum,INFO/ClippingRankSum,INFO/DB,INFO/DP,INFO/FS,INFO/MQ,INFO/MQRankSum,INFO/QD,INFO/ReadPosRankSum,INFO/SOR,INFO/ExcessHet \
| bcftools norm \
  --check-ref x \
  -m - \
  --fasta-ref /data/gusev/USERS/rlc47/TCGA/refs/GRCh37.fa \
| bcftools +fill-tags - -- -t AN,AC,AF \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
| bcftools view \
  --trim-alt-alleles \
  -Oz -o $WRKDIR/data/HMF.RAS_loci.$contig.vcf.gz \
  --include 'AC > 0'
tabix -p vcf -f $WRKDIR/data/HMF.RAS_loci.$contig.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/merge_sample_vcfs.$contig.sh
  rm $WRKDIR/LSF/logs/merge_sample_vcfs.$contig.*
  bsub \
    -q big -sla miket_sc \
    -n 4 -R 'rusage[mem=16000]' -J HMF_merge_sample_vcfs_$contig \
    -o $WRKDIR/LSF/logs/merge_sample_vcfs.$contig.log \
    -e $WRKDIR/LSF/logs/merge_sample_vcfs.$contig.err \
    $WRKDIR/LSF/scripts/merge_sample_vcfs.$contig.sh
done
# 2. Collect median GQ of autosomal homalt biallelic SNP GTs per sample
for contig in $( seq 1 22 ); do
  cat << EOF > $WRKDIR/LSF/scripts/get_median_homalt_snp_GQs.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools view \
  -m 2 -M 2 --types snps \
  $WRKDIR/data/HMF.RAS_loci.$contig.vcf.gz \
| bcftools query \
  -f '[%SAMPLE\t%GQ\n]' \
  -i 'GT="AA"' \
| gzip -c \
> $WRKDIR/data/sample_info/HMF.median_homalt_snp_gq.all_data.$contig.tsv.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/get_median_homalt_snp_GQs.$contig.sh
  rm $WRKDIR/LSF/logs/get_median_homalt_snp_GQs.$contig.*
  bsub \
    -q short -sla miket_sc -J HMF_get_median_homalt_snp_GQs_$contig \
    -o $WRKDIR/LSF/logs/get_median_homalt_snp_GQs.$contig.log \
    -e $WRKDIR/LSF/logs/get_median_homalt_snp_GQs.$contig.err \
    $WRKDIR/LSF/scripts/get_median_homalt_snp_GQs.$contig.sh
done
cat $WRKDIR/data/sample_info/HMF.median_homalt_snp_gq.all_data.*.tsv.gz \
> $WRKDIR/data/sample_info/HMF.median_homalt_snp_gq.all_data.tsv.gz
$CODEDIR/utils/calc_sample_medians.R \
  $WRKDIR/data/sample_info/HMF.median_homalt_snp_gq.all_data.tsv.gz \
| gzip -c > $WRKDIR/data/sample_info/HMF.median_homalt_gqs.tsv.gz
# 3. Fill reference GQs per sample with sample-specific median homalt GQ
for contig in $( seq 1 22 ) X Y; do
    cat << EOF > $WRKDIR/LSF/scripts/fill_missing_gqs.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
$CODEDIR/scripts/data_processing/fill_missing_vcf_format_values.py \
  $WRKDIR/data/HMF.RAS_loci.$contig.vcf.gz \
  $WRKDIR/data/sample_info/HMF.median_homalt_gqs.tsv.gz \
  $WRKDIR/data/HMF.RAS_loci.$contig.GQs_filled.vcf.gz
tabix -p vcf -f $WRKDIR/data/HMF.RAS_loci.$contig.GQs_filled.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/fill_missing_gqs.$contig.sh
  rm $WRKDIR/LSF/logs/fill_missing_gqs.$contig.*
  bsub \
    -q big -sla miket_sc \
    -n 4 -R 'rusage[mem=16000]' -J HMF_fill_missing_gqs_$contig \
    -o $WRKDIR/LSF/logs/fill_missing_gqs.$contig.log \
    -e $WRKDIR/LSF/logs/fill_missing_gqs.$contig.err \
    $WRKDIR/LSF/scripts/fill_missing_gqs.$contig.sh
done
# 4. Merge cohort-wide VCFs across all chromosomes
for contig in $( seq 1 22 ) X Y; do
  echo $WRKDIR/data/HMF.RAS_loci.$contig.vcf.gz
done > $WRKDIR/data/HMF.combined_vcf_shards.list
bcftools concat \
  --file-list $WRKDIR/data/HMF.combined_vcf_shards.list \
| bcftools annotate \
  --header-lines <( tabix -H $WRKDIR/refs/simple_hg19_header.vcf.gz | grep 'SOURCE\|SVTYPE' ) \
  -O z -o $WRKDIR/data/HMF.RAS_loci.vcf.gz
tabix -p vcf -f $WRKDIR/data/HMF.RAS_loci.vcf.gz


### Perform PCA & ancestry/kinship inference based on 1kG ref panel data
# Note: assumes ref panel of LD-independent common SNPs has already been produced
# and moved to google cloud. See harmonize_PCs_ancestries_PRS.sh for details.
# 1. Download per-sample VCFs subset in Terra
#    Note: This requires sample metadata manually downloaded from Terra
col_idxs=$( head -n1 $WRKDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
            | sed 's/\t/\n/g' | awk '{ if ($1 ~ "common_snps_vcf") print NR }' | paste -s -d, )
if ! [ -e $WRKDIR/data/1000G_pruned_common_snp_vcfs ]; then
  mkdir $WRKDIR/data/1000G_pruned_common_snp_vcfs
fi
sed '1d' $WRKDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
| cut -f$col_idxs | sed 's/\t/\n/g' \
| gsutil -m cp -I $WRKDIR/data/1000G_pruned_common_snp_vcfs/
# 2. Merge VCFs across all samples & apply AF/call rate filters
find $WRKDIR/data/1000G_pruned_common_snp_vcfs/ -name "*vcf.gz" \
> $WRKDIR/data/HMF.1000G_pruned_common_snp_vcfs.list
cat << EOF > $WRKDIR/LSF/scripts/HMF_merge_1000G_pruned_common_snp_vcfs.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools merge \
  --file-list $WRKDIR/data/HMF.1000G_pruned_common_snp_vcfs.list \
  --no-version \
  --missing-to-ref \
  --threads 4 \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
| bcftools view \
  --samples-file $WRKDIR/data/sample_info/HMF.ALL.samples.list \
| bcftools +fill-tags - -- -t AF,F_MISSING \
| bcftools norm \
  --atomize \
  --check-ref x \
  --fasta-ref /data/gusev/USERS/rlc47/TCGA/refs/GRCh37.fa \
  --threads 4 \
| bcftools view \
  -i 'INFO/AF > 0.01 & AF < 0.99 & INFO/F_MISSING < 0.01' \
  --type snps -m 2 -M 2 \
  -Oz -o $WRKDIR/data/HMF.1000G_pruned_common_snps.vcf.gz
tabix -p vcf -f $WRKDIR/data/HMF.1000G_pruned_common_snps.vcf.gz
EOF
chmod a+x $WRKDIR/LSF/scripts/HMF_merge_1000G_pruned_common_snp_vcfs.sh
rm $WRKDIR/LSF/logs/HMF_merge_1000G_pruned_common_snp_vcfs.*
bsub \
  -q big -sla miket_sc \
  -n 4 -R 'rusage[mem=16000]' -J HMF_merge_1000G_pruned_common_snp_vcfs \
  -o $WRKDIR/LSF/logs/HMF_merge_1000G_pruned_common_snp_vcfs.log \
  -e $WRKDIR/LSF/logs/HMF_merge_1000G_pruned_common_snp_vcfs.err \
  $WRKDIR/LSF/scripts/HMF_merge_1000G_pruned_common_snp_vcfs.sh
# 3. Further prune SNP markers with in-cohort frequencies from HMF
cat << EOF > $WRKDIR/LSF/scripts/prune_HMF_SNPs_inSampleFreq.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
module load plink/1.90b3
plink \
  --vcf $WRKDIR/data/HMF.1000G_pruned_common_snps.vcf.gz \
  --indep-pairwise 500kb 5 0.2 \
  --threads 4 \
  --memory 16000 \
  --out $WRKDIR/data/HMF.1000G_pruned_common_snps.inSample_pruned
bcftools view \
  -i "ID=@$WRKDIR/data/HMF.1000G_pruned_common_snps.inSample_pruned.prune.in" \
  -Oz -o $WRKDIR/data/HMF.1000G_pruned_common_snps.inSample_pruned.vcf.gz \
  $WRKDIR/data/HMF.1000G_pruned_common_snps.vcf.gz
tabix -p vcf -f $WRKDIR/data/HMF.1000G_pruned_common_snps.inSample_pruned.vcf.gz
EOF
chmod a+x $WRKDIR/LSF/scripts/prune_HMF_SNPs_inSampleFreq.sh
rm $WRKDIR/LSF/logs/prune_HMF_SNPs_inSampleFreq.*
bsub  -J prune_HMF_SNPs_inSampleFreq \
  -q big -sla miket_sc -n 4 -R 'rusage[mem=16000]' \
  -o $WRKDIR/LSF/logs/prune_HMF_SNPs_inSampleFreq.log \
  -e $WRKDIR/LSF/logs/prune_HMF_SNPs_inSampleFreq.err \
  $WRKDIR/LSF/scripts/prune_HMF_SNPs_inSampleFreq.sh
# 4. For ancestry truth labels, also extract genotypes for the same loci
# for all 1kG phase 3 samples and merge with HMF samples
for contig in $( seq 1 22 ); do
  echo $KGDIR/phase3/common_SNPs/1000G.ph3.common_snps.chr$contig.vcf.gz
done > $KGDIR/phase3/common_SNPs/1000G.ph3.common_snps.shards.list
cat << EOF > $WRKDIR/LSF/scripts/gather_1kG_SNPs_for_ancestry_inference.sh
bcftools concat \
  --file-list $KGDIR/phase3/common_SNPs/1000G.ph3.common_snps.shards.list \
  -a --regions-file $WRKDIR/data/HMF.1000G_pruned_common_snps.inSample_pruned.vcf.gz \
  --threads 4 \
| bcftools view \
  --no-version \
  -i "ID=@$WRKDIR/data/HMF.1000G_pruned_common_snps.inSample_pruned.prune.in" \
  -Oz -o $WRKDIR/data/1000G.ph3.pruned_common_snps.HMF_inSample_pruned.for_ancestry_inference.vcf.gz
tabix -p vcf -f $WRKDIR/data/1000G.ph3.pruned_common_snps.HMF_inSample_pruned.for_ancestry_inference.vcf.gz
bcftools merge \
  $WRKDIR/data/HMF.1000G_pruned_common_snps.inSample_pruned.vcf.gz \
  $WRKDIR/data/1000G.ph3.pruned_common_snps.HMF_inSample_pruned.for_ancestry_inference.vcf.gz \
  --threads 4 \
| bcftools +fill-tags - -- -t AF,F_MISSING,HWE \
| bcftools view \
  --include 'INFO/AF > 0.01 & INFO/AF < 0.99 & INFO/F_MISSING < 0.001 & INFO/HWE > 0.0000005' \
  -Oz -o $WRKDIR/data/HMF_plus_1kG.1000G_pruned_common_snps.inSample_pruned.vcf.gz
EOF
chmod a+x $WRKDIR/LSF/scripts/gather_1kG_SNPs_for_ancestry_inference.sh
rm $WRKDIR/LSF/logs/gather_1kG_SNPs_for_ancestry_inference.*
bsub  -J gather_1kG_SNPs_for_ancestry_inference \
  -q big -sla miket_sc -n 4 -R 'rusage[mem=16000]' \
  -o $WRKDIR/LSF/logs/gather_1kG_SNPs_for_ancestry_inference.log \
  -e $WRKDIR/LSF/logs/gather_1kG_SNPs_for_ancestry_inference.err \
  $WRKDIR/LSF/scripts/gather_1kG_SNPs_for_ancestry_inference.sh
# 6. Compute PCA & relatedness with plink
cat << EOF > $WRKDIR/LSF/scripts/HMF_plus_1000G_PCA.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
module load plink/2.0a2.3
plink2 \
  --threads 8 \
  --memory 16000 \
  --vcf $WRKDIR/data/HMF_plus_1kG.1000G_pruned_common_snps.inSample_pruned.vcf.gz \
  --pca \
  --make-king-table \
  --king-table-filter 0.06 \
  --out $WRKDIR/data/HMF_plus_1kG.1000G_pruned_common_snps.inSample_pruned
EOF
chmod a+x $WRKDIR/LSF/scripts/HMF_plus_1000G_PCA.sh
rm $WRKDIR/LSF/logs/HMF_plus_1000G_PCA.*
bsub \
  -q big-multi -sla miket_sc \
  -n 8 -R 'rusage[mem=16000]' -J HMF_plus_1000G_PCA \
  -o $WRKDIR/LSF/logs/HMF_plus_1000G_PCA.log \
  -e $WRKDIR/LSF/logs/HMF_plus_1000G_PCA.err \
  $WRKDIR/LSF/scripts/HMF_plus_1000G_PCA.sh
# 7. Assign ancestries using PedSV script
cd $CODEDIR/../ && \
git clone git@github.com:vanallenlab/ped_germline_SV.git && \
cd -
Rscript -e "install.packages(c('bedr', '$CODEDIR/../ped_germline_SV/src/PedSV_0.0.2.tar.gz') \
                             lib='~/R/x86_64-pc-linux-gnu-library/3.6', \
                             type='source', repos=NULL)"
# Note: 1000G truth labels must be uploaded to $WRKDIR/refs
if [ -e $WRKDIR/data/ancestry_inference ]; then
  rm -rf $WRKDIR/data/ancestry_inference
fi
mkdir $WRKDIR/data/ancestry_inference
$CODEDIR/../ped_germline_SV/gatksv_scripts/assign_ancestry.R \
  --PCs $WRKDIR/data/HMF_plus_1kG.1000G_pruned_common_snps.inSample_pruned.eigenvec \
  --training-labels $WRKDIR/refs/1000G_HGDP_training_labels.tsv.gz \
  --out-prefix $WRKDIR/data/ancestry_inference/HMF_ancestry_inference \
  --use-N-PCs 4 \
  --min-probability 0 \
  --plot
# 8. Update HMF sample metadata with new PCs and ancestry labels
$CODEDIR/scripts/data_processing/update_metadata_with_pcs_ancestry.py \
  --metadata $WRKDIR/data/sample_info/HMF.ALL.sample_metadata.tsv.gz \
  --pop-labels $WRKDIR/data/ancestry_inference/HMF_ancestry_inference.ancestry_labels.tsv \
  --PCs $WRKDIR/data/HMF_plus_1kG.1000G_pruned_common_snps.inSample_pruned.eigenvec

