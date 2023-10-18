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


### Curate germline variants in regions of interest
# 1. Copy single-sample VCFs from GCP to ERISOne
#    This requires sample metadata manually downloaded from Terra
col_idxs=$( head -n1 $WRKDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
            | sed 's/\t/\n/g' | awk '{ if ($1 ~ "ras_loci_vcf") print NR }' | paste -s -d, )
sed '1d' $WRKDIR/data/sample_info/HMF.hg19_terra_workspace.sample_info.tsv \
| cut -f$col_idxs | sed 's/\t/\n/g' \
| gsutil -m cp -I $WRKDIR/data/sample_vcfs/
# 2. Get consensus header from all samples
find $WRKDIR/data/sample_vcfs/ -name "*vcf.gz" > $WRKDIR/data/HMF.sample_vcfs.list
while read vcf; do
  tabix -H $vcf | grep '^##FILTER\|^##INFO'
done < $WRKDIR/data/HMF.sample_vcfs.list \
| sort -V | uniq
# 3. Merge single-sample VCFs into cohort-wide VCFs per chromosome
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
# 4. Merge cohort-wide VCFs across all chromosomes
for contig in $( seq 1 22 ) X Y; do
  echo $WRKDIR/data/HMF.RAS_loci.$contig.vcf.gz
done > $WRKDIR/data/HMF.combined_vcf_shards.list
bcftools concat \
  --file-list $WRKDIR/data/HMF.combined_vcf_shards.list \
  -O z -o $WRKDIR/data/HMF.RAS_loci.vcf.gz
tabix -p vcf -f $WRKDIR/data/HMF.RAS_loci.vcf.gz


### Curate somatic variants in regions of interest
# 1. Download & extract HMF somatic data
gsutil -m cp \
  gs://hmf-dr-355-us-central1/somatic.tar.gz \
  $WRKDIR/data/all_somatic/
tar -xzvf $WRKDIR/data/all_somatic/somatic.tar.gz





