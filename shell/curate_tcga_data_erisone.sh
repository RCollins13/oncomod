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
export CODEDIR=$WRKDIR/../code/ras_modifiers
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/sample_info LSF LSF/scripts LSF/logs refs; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


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
cut -f2 -d\  $GTDIR/TCGA.NORMAL.fam | sort | uniq \
> $WRKDIR/data/sample_info/TCGA.vcf.array_typed.samples.list
cut -f1-3 -d\- $WRKDIR/data/sample_info/TCGA.vcf.array_typed.samples.list | sort | uniq \
> $WRKDIR/data/sample_info/TCGA.vcf.array_typed.donors.list
# Get list of all IDs present in imputed VCF
tabix -H $GTDIR/IMPUTED/1.vcf.gz \
| fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
> $WRKDIR/data/sample_info/TCGA.vcf.array_imputed.samples.gusev_IDs.list
# Make map of Sasha's array IDs to canonical TCGA IDs
cut --complement -f2 -d\- $WRKDIR/data/sample_info/TCGA.vcf.array_imputed.samples.gusev_IDs.list \
> $WRKDIR/data/sample_info/TCGA.vcf.array_imputed.samples.list
cut -f1-3 -d\- $WRKDIR/data/sample_info/TCGA.vcf.array_imputed.samples.list | sort | uniq \
> $WRKDIR/data/sample_info/TCGA.vcf.array_imputed.donors.list
# Intersect samples with both exomes and genotyped arrays (not necessarily imputed)
$CODEDIR/scripts/data_processing/harmonize_tcga_samples.py \
  --exome-ids $WRKDIR/data/sample_info/TCGA.vcf.exome.samples.list \
  --exome-id-map $WRKDIR/refs/VanAllen.TCGA.WES_DeepVariant.sample_manifest.tsv.gz \
  --array-typed-ids $WRKDIR/data/sample_info/TCGA.vcf.array_typed.samples.list \
  --array-imputed-ids $WRKDIR/data/sample_info/TCGA.vcf.array_imputed.samples.gusev_IDs.list \
  --tcga-tss-table $CODEDIR/refs/TCGA_TSS_codes.tsv.gz \
  --tcga-study-table $CODEDIR/refs/TCGA_study_codes.tsv.gz \
  --out-prefix $WRKDIR/data/sample_info/TCGA


### Curate clinical information for patients of interest
# TODO: implement this
# $CODEDIR/scripts/data_processing/preprocess_dfci_profile_ehr.py \
#   --id-map-tsv $CLINDIR/TCGA_MRN_BL_PANEL.PBP.tab \
#   --genomic-csv $CLINDIR/OncDRS/ALL_2021_11/GENOMIC_SPECIMEN.csv.gz \
#   --dx-csv $CLINDIR/OncDRS/ALL_2021_11/CANCER_DIAGNOSIS_CAREG.csv.gz \
#   --ancestry-csv $CLINDIR/TCGA_2022_ANCESTRY.csv.gz \
#   --hx-csv $CLINDIR/OncDRS/ALL_2021_11/HEALTH_HISTORY.csv.gz \
#   --survival-csv $CLINDIR/OncDRS/ALL_2021_11/PT_INFO_STATUS_REGISTRATION.csv.gz \
#   --out-prefix $WRKDIR/data/sample_info/TCGA. \
#   --vcf-ids $WRKDIR/data/sample_info/TCGA.vcf.samples.list


### Subset VCFs to patients of interest and RAS loci
# Extract samples & loci of interest from Sasha's imputed arrays
while read contig start end gene; do
  cat << EOF > $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools view \
  -O z -o $WRKDIR/data/TCGA.$gene.array.vcf.gz \
  --min-ac 1 \
  --samples-file $WRKDIR/data/sample_info/TCGA.ALL.array.samples.list \
  --regions "$contig:${start}-$end" \
  $GTDIR/$contig.vcf.gz
tabix -p vcf -f $WRKDIR/data/TCGA.$gene.array.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
  rm $WRKDIR/LSF/logs/extract_${gene}_variants.*
  bsub \
    -q long -R 'rusage[mem=6000]' -n 2 -J TCGA_extract_${gene}_variants \
    -o $WRKDIR/LSF/logs/extract_${gene}_variants.log \
    -e $WRKDIR/LSF/logs/extract_${gene}_variants.err \
    $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" )
# Define well-covered exome intervals within RAS loci of interest
$TMPDIR/define_well_covered_targets.py \
  --coverage-matrix $WRKDIR/data/101122_TCGA_10960_interval_process.tsv.gz \
  --samples-list $WRKDIR/data/sample_info/TCGA.ALL.exome.samples.list \
  --min-frac-samples 0.9 \
  --min-frac-target 0.9
# TODO: intersect well-covered intervals with RAS gene loci to get regions for extraction
# Use bcftools to stream WES data from gs:// bucket for samples & loci of interest
gsutil -m cp \
  gs://fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz.tbi \
  $WRKDIR/data/
gcloud auth application-default login
export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
cd $WRKDIR/data
while read contig start end gene; do
  bcftools view \
    -O z -o $WRKDIR/data/TCGA.$gene.exome.vcf.gz \
    --min-ac 1 \
    --samples-file $WRKDIR/data/sample_info/TCGA.ALL.exome.samples.list \
    --regions "$contig:${start}-$end" \
    gs://fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz
  tabix -p vcf -f $WRKDIR/data/TCGA.$gene.exome.vcf.gz
done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" )
# Merge VCFs for exomes and arrays for each gene
# TODO: implement this
# Merge VCFs for eacn gene into a single VCF and index the merged VCF
zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 \
| xargs -I {} echo "$WRKDIR/data/TCGA.{}.vcf.gz" \
> $WRKDIR/data/TCGA.single_gene_vcfs.list
bcftools concat \
  --file-list $WRKDIR/data/TCGA.single_gene_vcfs.list \
  -O z -o $WRKDIR/data/TCGA.RAS_loci.vcf.gz
tabix -p vcf -f $WRKDIR/data/TCGA.RAS_loci.vcf.gz

