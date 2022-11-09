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
for SUBDIR in data data/sample_info data/sample_info/TCGA_BMI LSF LSF/scripts LSF/logs refs misc; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Recode array-genotyped SNPs as VCF
module load plink/1.90b3
plink \
  --bfile $GTDIR/TCGA.NORMAL \
  --recode vcf bgz \
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
  --tcga-tss-table $CODEDIR/refs/TCGA_TSS_codes.tsv.gz \
  --tcga-study-table $CODEDIR/refs/TCGA_study_codes.tsv.gz \
  --out-prefix $WRKDIR/data/sample_info/TCGA


### Curate clinical information for patients of interest
# Download & extract TCGA BMI
if [ -e $WRKDIR/data/sample_info/TCGA.BMI.tsv ]; then
  rm $WRKDIR/data/sample_info/TCGA.BMI.tsv
fi
for prefix in PAAD COAD READ SKCM; do
  wget \
    -P $WRKDIR/data/sample_info/TCGA_BMI/ \
    https://hgdownload.soe.ucsc.edu/gbdb/hg38/gdcCancer/$prefix.bb
  /data/talkowski/tools/bin/bigBedToBed \
    $WRKDIR/data/sample_info/TCGA_BMI/$prefix.bb \
    /dev/stdout \
  | awk -v FS="\t" -v OFS="\t" '{ if ($29!="--") print $35, $29 }' \
  | $CODEDIR/scripts/data_processing/parse_ucsc_tcga_bmi.py \
  | sort -Vk1,1 | uniq \
  >> $WRKDIR/data/sample_info/TCGA_BMI/TCGA.BMI.tsv
  rm $WRKDIR/data/sample_info/TCGA_BMI/$prefix.bb
done
gzip -f $WRKDIR/data/sample_info/TCGA_BMI/TCGA.BMI.tsv
# Note: TCGA ancestry label assignments came from Carrot-Zhang et al., Cancer Cell, 2020
# https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020
# Filename: Broad_ancestry_PCA.txt
# This file has been renamed and relocated to $WRKDIR/data/TCGA.ancestry.tsv.gz
$CODEDIR/scripts/data_processing/preprocess_tcga_phenotypes.py \
  --id-map-tsv $WRKDIR/data/sample_info/TCGA.ALL.id_map.tsv.gz \
  --cdr-csv $BASEDIR/TCGA_CDR.csv \
  --tcga-study-table $CODEDIR/refs/TCGA_study_codes.tsv.gz \
  --bmi-tsv $WRKDIR/data/sample_info/TCGA_BMI/TCGA.BMI.tsv.gz \
  --ancestry-tsv $WRKDIR/data/TCGA.ancestry.tsv.gz \
  --pcs-txt $GTDIR/TCGA.COMBINED.QC.NORMAL.eigenvec \
  --out-prefix $WRKDIR/data/sample_info/TCGA.


### Subset VCFs to patients of interest and RAS loci
# Extract samples & loci of interest from genotyped and imputed arrays
while read contig start end gene; do
  # for tech in array_typed array_imputed; do
  for tech in array_imputed; do
    bcftools_options="--min-ac 1"
    case $tech in
      array_typed)
        VCF=$WRKDIR/data/TCGA.array_typed.vcf.gz
        ;;
      array_imputed)
        VCF=$GTDIR/IMPUTED/$contig.vcf.gz
        bcftools_options="$bcftools_options --exclude 'INFO/INFO < 0.8'"
        ;;
    esac
    cat << EOF > $WRKDIR/LSF/scripts/extract_${gene}_variants_${tech}.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools view \
  -O z -o $WRKDIR/data/TCGA.$gene.$tech.vcf.gz \
  $bcftools_options \
  --samples-file $WRKDIR/data/sample_info/TCGA.ALL.$tech.samples.list \
  --regions "$contig:${start}-$end" \
  $VCF
tabix -p vcf -f $WRKDIR/data/TCGA.$gene.$tech.vcf.gz
EOF
    chmod a+x $WRKDIR/LSF/scripts/extract_${gene}_variants_${tech}.sh
    rm $WRKDIR/LSF/logs/extract_${gene}_variants_${tech}.*
    bsub \
      -q short -R 'rusage[mem=6000]' -n 2 -J TCGA_extract_${gene}_variants_${tech} \
      -o $WRKDIR/LSF/logs/extract_${gene}_variants_${tech}.log \
      -e $WRKDIR/LSF/logs/extract_${gene}_variants_${tech}.err \
      $WRKDIR/LSF/scripts/extract_${gene}_variants_${tech}.sh
  done
done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" )
# Define well-covered exome intervals within RAS loci of interest
$CODEDIR/scripts/data_processing/define_well_covered_targets.py \
  --coverage-matrix $WRKDIR/data/101122_TCGA_10960_interval_process.tsv.gz \
  --samples-list $WRKDIR/data/sample_info/TCGA.ALL.exome.samples.list \
  --min-frac-samples 0.9 \
  --min-frac-target 0.9 \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> $WRKDIR/refs/TCGA_WES.covered_intervals.bed.gz
tabix -p bed -f $WRKDIR/refs/TCGA_WES.covered_intervals.bed.gz
while read contig start end gene; do
  bedtools intersect -u \
    -a $WRKDIR/refs/TCGA_WES.covered_intervals.bed.gz \
    -b <( echo -e "$contig\t$start\t$end" ) \
  | bgzip -c > $WRKDIR/refs/TCGA_WES.covered_intervals.$gene.bed.gz
  tabix -p bed -f $WRKDIR/refs/TCGA_WES.covered_intervals.$gene.bed.gz
done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" )
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
    --regions-file $WRKDIR/refs/TCGA_WES.covered_intervals.$gene.bed.gz \
    gs://fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz
  tabix -p vcf -f $WRKDIR/data/TCGA.$gene.exome.vcf.gz
done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" )
# Merge VCFs for exomes and arrays for each gene
while read gene; do
    cat << EOF > $WRKDIR/LSF/scripts/merge_arrays_exomes_${gene}.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
$CODEDIR/scripts/data_processing/merge_tcga_arrays_exomes.py \
  --sample-id-map /data/gusev/USERS/rlc47/TCGA/data/sample_info/TCGA.ALL.id_map.tsv.gz \
  --exome-vcf $WRKDIR/data/TCGA.$gene.exome.vcf.gz \
  --array-typed-vcf $WRKDIR/data/TCGA.$gene.array_typed.vcf.gz \
  --array-imputed-vcf $WRKDIR/data/TCGA.$gene.array_imputed.vcf.gz \
  --ref-fasta $WRKDIR/refs/GRCh37.fa \
  --header $WRKDIR/refs/simple_hg19_header.vcf.gz \
  --outfile $WRKDIR/data/TCGA.$gene.merged.vcf.gz \
  --verbose
tabix -p vcf -f $WRKDIR/data/TCGA.$gene.merged.vcf.gz
EOF
    chmod a+x $WRKDIR/LSF/scripts/merge_arrays_exomes_${gene}.sh
    rm $WRKDIR/LSF/logs/merge_arrays_exomes_${gene}.*
    bsub \
      -q normal -R 'rusage[mem=6000]' -J TCGA_merge_arrays_exomes_${gene} \
      -o $WRKDIR/LSF/logs/merge_arrays_exomes_${gene}.log \
      -e $WRKDIR/LSF/logs/merge_arrays_exomes_${gene}.err \
      $WRKDIR/LSF/scripts/merge_arrays_exomes_${gene}.sh
done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 )
# Merge VCFs for eacn gene into a single VCF and index the merged VCF
zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 \
| xargs -I {} echo "$WRKDIR/data/TCGA.{}.merged.vcf.gz" \
> $WRKDIR/data/TCGA.single_gene_vcfs.list
bcftools concat \
  --file-list $WRKDIR/data/TCGA.single_gene_vcfs.list \
  -O z -o $WRKDIR/data/TCGA.RAS_loci.vcf.gz
tabix -p vcf -f $WRKDIR/data/TCGA.RAS_loci.vcf.gz


### Curate somatic data for patients of interest
# Preprocess mc3 data by subsetting on samples of interest to improve pandas IO speed
cat \
  <( zcat $WRKDIR/data/mc3.v0.2.8.PUBLIC.maf.gz | head -n1 ) \
  <( zcat $WRKDIR/data/mc3.v0.2.8.PUBLIC.maf.gz \
     | fgrep -wf $WRKDIR/data/sample_info/TCGA.ALL.donors.list ) \
| gzip -c > $WRKDIR/data/mc3.v0.2.8.PUBLIC.SKCM_PDAC_CRAD.maf.gz
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
zcat $WRKDIR/data/mc3.v0.2.8.PUBLIC.SKCM_PDAC_CRAD.maf.gz \
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
  --mc3-tsv $WRKDIR/data/mc3.v0.2.8.PUBLIC.SKCM_PDAC_CRAD.maf.gz \
  --cna-bed $WRKDIR/data/TCGA.CNA.b19.bed.gz \
  --donors-list $WRKDIR/data/sample_info/TCGA.ALL.donors.list \
  --purity-tsv $WRKDIR/data/TCGA.tumor_purity.Aran_2015.tsv.gz \
  --genes-gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --priority-genes $WRKDIR/../refs/COSMIC.all_GCG.Nov8_2022.genes.list \
  --outfile $WRKDIR/data/TCGA.somatic_variants.tsv.gz


# Summarize somatic variant status by gene & cancer type
for cancer in PDAC CRAD SKCM; do
  while read gene; do
    n_samp=$( fgrep -wvf \
                $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.list \
                ${WRKDIR}/data/sample_info/TCGA.$cancer.donors.list | wc -l )
    zcat $WRKDIR/data/TCGA.somatic_variants.tsv.gz \
    | fgrep -wf ${WRKDIR}/data/sample_info/TCGA.$cancer.donors.list \
    | awk -v gene=$gene -v FS="\t" '{ if ($8==gene && ($14~/AMP|DEL|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/)) print $6 }' \
    | sort | uniq | wc -l | awk -v n=$n_samp '{ print $1/n }'
  done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 ) | paste -s -
done

