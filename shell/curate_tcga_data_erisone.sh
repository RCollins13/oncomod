#!/usr/bin/env bash

################################
#    EGFR Modifiers Project    #
################################

# Copyright (c) 2023-Present Ryan L. Collins, Jackie LoPiccolo, and the Gusev/Van Allen Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Curate (some) data from TCGA cohort for RAS modifiers study

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export BASEDIR=/data/gusev/TCGA
export GTDIR=/data/gusev/TCGA/GENOTYPES
export WRKDIR=/data/gusev/USERS/rlc47/TCGA
export CODEDIR=$WRKDIR/../EGFR_modifier_analysis/code/ras_modifiers
cd $WRKDIR


### NOTE: EXPECTS FULL RAS MODIFIER DATA CURATION TO BE ALREADY COMPLETE 
### SEE VERSION OF THIS FILE IN THE MAIN BRANCH


### Subset VCFs to patients of interest and RAS loci
# Extract samples & loci of interest from genotyped and imputed arrays
while read contig start end gene; do
  for tech in array_typed array_imputed; do
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
done < <( zcat $CODEDIR/refs/EGFR_loci.GRCh37.bed.gz | fgrep -v "#" )
# Define well-covered exome intervals within EGFR locus of interest
while read contig start end gene; do
  bedtools intersect -u \
    -a $WRKDIR/refs/TCGA_WES.covered_intervals.bed.gz \
    -b <( echo -e "$contig\t$start\t$end" ) \
  | bgzip -c > $WRKDIR/refs/TCGA_WES.covered_intervals.$gene.bed.gz
  tabix -p bed -f $WRKDIR/refs/TCGA_WES.covered_intervals.$gene.bed.gz
done < <( zcat $CODEDIR/refs/EGFR_loci.GRCh37.bed.gz | fgrep -v "#" )
# Use bcftools to stream WES data from gs:// bucket for samples & loci of interest
gsutil -m cp \
  gs://terra-workspace-archive-lifecycle/fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz.tbi \
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
    gs://terra-workspace-archive-lifecycle/fc-e5ae96e4-c495-44d1-9155-b27057d570d8/e3d12a6b-5051-4656-b911-ea425aa14ce7/VT_Decomp/49610c65-a66c-45e9-9ef3-a3b5ab80ac9f/call-VTRecal/all_normal_samples.vt2_normalized_spanning_alleles.vcf.gz
  tabix -p vcf -f $WRKDIR/data/TCGA.$gene.exome.vcf.gz
done < <( zcat $CODEDIR/refs/EGFR_loci.GRCh37.bed.gz | fgrep -v "#" )
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
done < <( zcat $CODEDIR/refs/EGFR_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 )
# Merge VCFs for eacn gene into a single VCF and index the merged VCF
zcat $CODEDIR/refs/EGFR_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 \
| xargs -I {} echo "$WRKDIR/data/TCGA.{}.merged.vcf.gz" \
> $WRKDIR/data/TCGA.single_gene_vcfs.list
bcftools concat \
  --file-list $WRKDIR/data/TCGA.single_gene_vcfs.list \
  -O z -o $WRKDIR/data/TCGA.EGFR_loci.vcf.gz
tabix -p vcf -f $WRKDIR/data/TCGA.EGFR_loci.vcf.gz


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
  --priority-genes $WRKDIR/../refs/COSMIC.all_GCG.Nov8_2022.genes.list \
  --ref-fasta $WRKDIR/refs/GRCh37.fa \
  --header $WRKDIR/../refs/simple_hg19_header.somatic.vcf.gz \
  --outfile $WRKDIR/data/TCGA.somatic_variants.vcf.gz \
  --out-tsv $WRKDIR/data/TCGA.somatic_variants.tsv.gz
tabix -p vcf -f $WRKDIR/data/TCGA.somatic_variants.vcf.gz


# Summarize somatic variant status by gene & cancer type
for cancer in PDAC CRAD LUAD SKCM; do
  n_samp=$( fgrep -wvf \
              $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.list \
              ${WRKDIR}/data/sample_info/TCGA.$cancer.donors.list | wc -l )
  while read gene; do
    zcat $WRKDIR/data/TCGA.somatic_variants.tsv.gz \
    | fgrep -wf ${WRKDIR}/data/sample_info/TCGA.$cancer.donors.list \
    | fgrep -wvf $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.list \
    | awk -v gene=$gene -v FS="\t" '{ if ($9==gene && ($15~/AMP|DEL|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/)) print $7 }' \
    | sort | uniq | wc -l | awk -v n=$n_samp '{ print $1/n }'
  done < <( zcat $CODEDIR/refs/EGFR_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 )
  zcat $WRKDIR/data/TCGA.somatic_variants.tsv.gz \
    | fgrep -wf ${WRKDIR}/data/sample_info/TCGA.$cancer.donors.list \
    | fgrep -wvf $WRKDIR/data/sample_info/TCGA.ALL.donors.missing_somatic.list \
    | awk -v FS="\t" '{ if ($9~/NRAS|HRAS|KRAS/ && $15~/AMP|DEL|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/) print $7 }' \
    | sort | uniq | wc -l | awk -v n=$n_samp '{ print $1/n }'
done | paste - - - -

