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
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/sample_vcfs data/sample_info data/all_somatic LSF LSF/scripts LSF/logs ../refs refs; do
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
# $CODEDIR/scripts/data_processing/preprocess_hmf_phenotypes.py \
$TMPDIR/preprocess_hmf_phenotypes.py \
  --metadata $WRKDIR/data/sample_info/metadata.tsv \
  --purple-qc $WRKDIR/data/all_somatic/HMF.all.purple.somatic_qc.tsv.gz \
  --vcf-ids $WRKDIR/data/HMF.ss_germline_vcfs.all_samples.list \
  --hmf-json $WRKDIR/data/sample_info/manifest.json \
  --out-prefix $WRKDIR/data/sample_info/HMF.


### Curate somatic variants in regions of interest
# 2. Clear data for unnecessary samples
# TODO: implement this



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
| bcftools +fill-tags - -- -t AN,AC,AF,HWE \
| bcftools view \
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
# 2. Merge cohort-wide VCFs across all chromosomes
for contig in $( seq 1 22 ) X Y; do
  echo $WRKDIR/data/HMF.RAS_loci.$contig.vcf.gz
done > $WRKDIR/data/HMF.combined_vcf_shards.list
bcftools concat \
  --file-list $WRKDIR/data/HMF.combined_vcf_shards.list \
  -O z -o $WRKDIR/data/HMF.RAS_loci.vcf.gz
tabix -p vcf -f $WRKDIR/data/HMF.RAS_loci.vcf.gz






