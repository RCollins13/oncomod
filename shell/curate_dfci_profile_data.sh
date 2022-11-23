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
export CODEDIR=$WRKDIR/../code/ras_modifiers
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/sample_info LSF LSF/scripts LSF/logs ../refs refs; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Curate clinical information for patients of interest
# Get list of all IDs present in VCF
tabix -H $GTDIR/PROFILE_COMB.22.HQ.vcf.gz \
| fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
> $WRKDIR/data/sample_info/PROFILE.vcf.samples.list
# Curate EHR and get list of patients from cancer types of interest
$CODEDIR/scripts/data_processing/preprocess_dfci_profile_ehr.py \
  --id-map-tsv $CLINDIR/PROFILE_MRN_BL_PANEL.PBP.tab \
  --genomic-csv $CLINDIR/OncDRS/ALL_2021_11/GENOMIC_SPECIMEN.csv.gz \
  --dx-csv $CLINDIR/OncDRS/ALL_2021_11/CANCER_DIAGNOSIS_CAREG.csv.gz \
  --ancestry-csv $CLINDIR/PROFILE_2022_ANCESTRY.csv.gz \
  --hx-csv $CLINDIR/OncDRS/ALL_2021_11/HEALTH_HISTORY.csv.gz \
  --survival-csv $CLINDIR/OncDRS/ALL_2021_11/PT_INFO_STATUS_REGISTRATION.csv.gz \
  --out-prefix $WRKDIR/data/sample_info/PROFILE. \
  --vcf-ids $WRKDIR/data/sample_info/PROFILE.vcf.samples.list


### Subset VCFs to patients of interest and RAS loci
# Recompress chromosomes 1-3 with bgzip & reindex
for contig in $( seq 1 3 ); do
  bsub -R 'rusage[mem=8000]' -q long -J compress_index_$contig \
    "cd /data/gusev/PROFILE/2020_2022_combined/IMPUTE_HQ; \
     . /PHShome/rlc47/.bashrc; \
     zcat PROFILE_COMB.$contig.HQ.vcf.gz | bgzip -c > PROFILE_COMB.$contig.HQ.vcf.bgz; \
     bcftools index PROFILE_COMB.$contig.HQ.vcf.bgz"
done
# Extract samples & loci of interest
while read contig start end gene; do
  cat << EOF > $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools view \
  -O z -o $WRKDIR/data/PROFILE.$gene.vcf.gz \
  --min-ac 1 \
  --samples-file $WRKDIR/data/sample_info/PROFILE.ALL.samples.list \
  --regions "$contig:${start}-$end" \
  $GTDIR/PROFILE_COMB.$contig.HQ.vcf.gz
tabix -p vcf -f $WRKDIR/data/PROFILE.$gene.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
  rm $WRKDIR/LSF/logs/extract_${gene}_variants.*
  bsub \
    -q long -R 'rusage[mem=6000]' -n 2 -J PROFILE_extract_${gene}_variants \
    -o $WRKDIR/LSF/logs/extract_${gene}_variants.log \
    -e $WRKDIR/LSF/logs/extract_${gene}_variants.err \
    $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" )
# Merge VCFs for each gene into a single VCF and index the merged VCF
zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 \
| xargs -I {} echo "$WRKDIR/data/PROFILE.{}.vcf.gz" \
> $WRKDIR/data/PROFILE.single_gene_vcfs.list
bcftools concat \
  --file-list $WRKDIR/data/PROFILE.single_gene_vcfs.list \
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
$TMPDIR/preprocess_dfci_profile_somatic.py \
  --mutation-csv $CLINDIR/OncDRS/ALL_2021_11/GENOMIC_MUTATION_RESULTS.csv.gz \
  --cna-csv $CLINDIR/OncDRS/ALL_2021_11/GENOMIC_CNV_RESULTS.csv.gz \
  --samples-list $WRKDIR/data/sample_info/PROFILE.ALL.samples.list \
  --no-mutation-data $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic.mut.list \
  --no-cna-data $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic.cna.list \
  --id-map-tsv $CLINDIR/PROFILE_MRN_BL_PANEL.PBP.tab \
  --genes-gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --ref-fasta $WRKDIR/refs/GRCh37.fa \
  --header $WRKDIR/../refs/simple_hg19_header.somatic.vcf.gz \
  --outfile $WRKDIR/data/PROFILE.somatic_variants.vcf.gz


# Summarize somatic variant status by gene & cancer type
for cancer in PDAC CRAD SKCM; do
  n_samp=$( fgrep -wvf \
              $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic.list \
              ${WRKDIR}/data/sample_info/PROFILE.$cancer.samples.list | wc -l )
  while read gene; do
    zcat $WRKDIR/data/PROFILE.somatic_variants.tsv.gz \
    | fgrep -wf ${WRKDIR}/data/sample_info/PROFILE.$cancer.samples.list \
    | fgrep -wvf $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic.list \
    | awk -v gene=$gene -v FS="\t" '{ if ($7==gene && ($14~/AMP|DEL|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/)) print $6 }' \
    | sort | uniq | wc -l | awk -v n=$n_samp '{ print $1/n }'
  done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 )
  zcat $WRKDIR/data/PROFILE.somatic_variants.tsv.gz \
  | fgrep -wf ${WRKDIR}/data/sample_info/PROFILE.$cancer.samples.list \
  | fgrep -wvf $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic.list \
  | awk -v FS="\t" '{ if ($7~/NRAS|HRAS|KRAS/ && ($14~/AMP|DEL|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/)) print $6 }' \
  | sort | uniq | wc -l | awk -v n=$n_samp '{ print $1/n }'
done | paste - - - -

