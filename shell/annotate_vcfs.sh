#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Annotate all VCFs prior to analysis

# Note: intended to be executed on the MGB ERISOne cluster


##################
### Basic prep ###
##################

### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export CODEDIR=$WRKDIR/../code/oncomod
export VEP_CACHE=$WRKDIR/../refs/vep_cache
export VEP_PLUGINS=$WRKDIR/../code/vep_plugins
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in LSF LSF/logs LSF/scripts; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done
if ! [ -e $VEP_CACHE ]; then
  mkdir $VEP_CACHE
fi
if ! [ -e $VEP_PLUGINS ]; then
  mkdir $VEP_PLUGINS
fi


### Download & install VEP and associated plugins
# Prep directory structure
cd $WRKDIR/../code/ && \
git clone https://github.com/Ensembl/ensembl-vep.git && \
cd ensembl-vep && \
perl INSTALL.pl \
  --NO_HTSLIB \
  --CONVERT \
  --CACHEDIR $VEP_CACHE/ \
  --PLUGINSDIR $VEP_PLUGINS/ && \
cd $WRKDIR
# Download GRCh37 & GRCh38 indexed VEP caches
for ref in GRCh37 GRCh38; do
  bsub \
    -q normal -J ${ref}_VEP_cache_download \
    -o $WRKDIR/LSF/logs/${ref}_VEP_cache_download.log \
    -e $WRKDIR/LSF/logs/${ref}_VEP_cache_download.err \
    "cd $VEP_CACHE && \
     wget https://ftp.ensembl.org/pub/release-108/variation/indexed_vep_cache/homo_sapiens_vep_108_${ref}.tar.gz && \
     tar -xzf homo_sapiens_vep_108_${ref}.tar.gz"
done
# Test VEP installation
cd $WRKDIR/../code/ensembl-vep && \
./vep \
  -i examples/homo_sapiens_GRCh37.vcf \
  --cache \
  --dir_cache $VEP_CACHE/ \
  --cache_version 108 \
  --offline && \
cd $WRKDIR


#####################
### Data curation ###
#####################

### Curate custom datasets for use by VEP
# Prep subdirectories
for subdir in UTRAnnotator gnomad loftee CADD clinvar cosmic SpliceAI; do
  if ! [ -e $VEP_CACHE/$subdir ]; then
    mkdir $VEP_CACHE/$subdir
  fi
done
# dbNSFP
bsub \
  -q normal -J dbNSFP_cache_download \
  -o $WRKDIR/LSF/logs/dbNSFP_cache_download.log \
  -e $WRKDIR/LSF/logs/dbNSFP_cache_download.err \
  "cd $VEP_CACHE && \
   wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.3a.zip && \
   unzip dbNSFP4.3a.zip"
zcat $VEP_CACHE/dbNSFP4.3a_variant.chr1.gz \
| head -n1 > $TMPDIR/dbNSFP_header
cat <<EOF > $WRKDIR/LSF/scripts/prep_dbNSFP.grch38.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
zgrep -h -v ^#chr $VEP_CACHE/dbNSFP4.3a_variant.chr* \
| sort -T $TMPDIR -k1,1 -k2,2n - | cat $TMPDIR/dbNSFP_header - \
| bgzip -c > $VEP_CACHE/dbNSFP4.3a_grch38.gz
tabix -s 1 -b 2 -e 2 -f $VEP_CACHE/dbNSFP4.3a_grch38.gz
EOF
cat <<EOF > $WRKDIR/LSF/scripts/prep_dbNSFP.grch37.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
zgrep -h -v ^#chr $VEP_CACHE/dbNSFP4.3a_variant.chr* \
| awk '\$8 != "." ' | sort -T $TMPDIR -k8,8 -k9,9n - | cat $TMPDIR/dbNSFP_header - \
| bgzip -c > $VEP_CACHE/dbNSFP4.3a_grch37.gz
tabix -s 8 -b 9 -e 9 -f $VEP_CACHE/dbNSFP4.3a_grch37.gz
EOF
for ref in 37 38; do
  script=$WRKDIR/LSF/scripts/prep_dbNSFP.grch$ref.sh
  chmod a+x $script
  rm $WRKDIR/LSF/logs/prep_dbNSFP.grch$ref.*
  bsub \
    -q normal -J prep_dbNSFP.grch$ref \
    -o $WRKDIR/LSF/logs/prep_dbNSFP.grch$ref.log \
    -e $WRKDIR/LSF/logs/prep_dbNSFP.grch$ref.err \
    $WRKDIR/LSF/scripts/prep_dbNSFP.grch$ref.sh
done
# LOFTEE
for file in $VEP_PLUGINS/loftee/*; do
  ln -fs $( echo -e $file ) $VEP_PLUGINS/$( basename $file )
done
for suffix in gz gz.fai gz.gzi; do
  bsub \
    -q normal -J dl_loftee_ancestor_fa$suffix \
    -o $WRKDIR/LSF/logs/dl_loftee_ancestor_fa$suffix.log \
    -e $WRKDIR/LSF/logs/dl_loftee_ancestor_fa$suffix.err \
    "cd $VEP_CACHE/loftee && \
     wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.$suffix"
done
bsub \
  -q normal -J dl_loftee_conservation_db \
  -o $WRKDIR/LSF/logs/dl_loftee_conservation_db.log \
  -e $WRKDIR/LSF/logs/dl_loftee_conservation_db.err \
  "cd $VEP_CACHE/loftee && \
   wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz && \
   gunzip phylocsf_gerp.sql.gz"
# CADD
for prefix in whole_genome_SNVs InDels; do
  for suffix in gz gz.tbi; do
    bsub \
      -q normal -J dl_cadd_$prefix.tsv.$suffix \
      -o $WRKDIR/LSF/logs/dl_cadd_$prefix.tsv.$suffix.log \
      -e $WRKDIR/LSF/logs/dl_cadd_$prefix.tsv.$suffix.err \
      "cd $VEP_CACHE/CADD && \
       wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/$prefix.tsv.$suffix"
  done
done
# ClinVar
wget \
  -P $VEP_CACHE/clinvar/ \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
# COSMIC
# Note: first need to generate a COSMIC command line access key
# See https://cancer.sanger.ac.uk/cosmic/download for instructions
# Need to assign this to variable $my_access_code for these commands to work
curl \
  -H "Authorization: Basic $my_access_code" \
  https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/CMC.tar \
| jq .url | tr -d '"' \
| xargs -I {} wget --no-check-certificate -P $TMPDIR/ {}
mkdir $TMPDIR/cosmic/ && tar -xvf $TMPDIR/CMC.tar* -C $TMPDIR/cosmic/
bsub \
  -q normal -J cmc_to_vcf \
  -o $WRKDIR/LSF/logs/cmc_to_vcf.log \
  -e $WRKDIR/LSF/logs/cmc_to_vcf.err \
  "$CODEDIR/scripts/data_processing/cmc_to_vcf.py \
     --cmc $TMPDIR/cosmic/cmc_export.tsv.gz \
     --header $WRKDIR/../refs/simple_hg19_header.somatic.vcf.gz \
     --ref-fasta $TCGADIR/refs/GRCh37.fa \
   | bcftools sort -O z -o $VEP_CACHE/cosmic/cmc.vcf.gz && \
   tabix -p vcf -f $VEP_CACHE/cosmic/cmc.vcf.gz"
# ABC
wget -P $TMPDIR/ \
  ftp://ftp.broadinstitute.org/outgoing/lincRNA/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz
$CODEDIR/scripts/data_processing/prep_ABC.py \
  --abc $TMPDIR/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> $VEP_CACHE/ABC_enhancers.bed.gz
tabix -p bed -f $VEP_CACHE/ABC_enhancers.bed.gz
# gnomAD
for subset in exomes genomes; do
  bsub \
    -q normal -J gnomad_${subset}_cache_download \
    -o $WRKDIR/LSF/logs/gnomad_${subset}_cache_download.log \
    -e $WRKDIR/LSF/logs/gnomad_${subset}_cache_download.err \
    "cd $VEP_CACHE/gnomad && \
     wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.$subset.r2.0.1.sites.noVEP.vcf.gz && \
     wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.$subset.r2.0.1.sites.noVEP.vcf.gz.tbi"
done
# Gene promoters (defined as +2kb upstream of TSS)
$CODEDIR/scripts/data_processing/make_promoters.py \
  --gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --coding-only \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> $VEP_CACHE/gencode_promoters.bed.gz
tabix -p bed -f $VEP_CACHE/gencode_promoters.bed.gz
# CCRs
wget -P $TMPDIR/ \
  https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz \
  https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz
zcat \
  $TMPDIR/ccrs.autosomes.v2.20180420.bed.gz \
  $TMPDIR/ccrs.xchrom.v2.20180420.bed.gz \
| awk -v OFS="\t" '{ print $1, $2, $3, $5":"$4 }' \
| grep -ve '^#' | sort -Vk1,1 -k2,2n -k3,3 | bgzip -c \
> $VEP_CACHE/CCR.bed.gz
tabix -f -p bed $VEP_CACHE/CCR.bed.gz
# GERP, phastCons, phyloP
cat << EOF > $TMPDIR/conservation_bw_urls.list
http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw
EOF
while read url; do
  bsub \
    -q normal -J dl_$( basename $url ) \
    -o $WRKDIR/LSF/logs/dl_$( basename $url ).log \
    -e $WRKDIR/LSF/logs/dl_$( basename $url ).err \
    "wget -P $VEP_CACHE/ $url"
done < $TMPDIR/conservation_bw_urls.list
while read url; do
  old=$VEP_CACHE/$( basename $url )
  new=$( echo $old | sed 's/hg19/GRCh37/g' )
  bsub \
    -q normal -J convert_$( basename $old )_contigs \
    -o $WRKDIR/LSF/logs/convert_$( basename $old )_contigs.log \
    -e $WRKDIR/LSF/logs/convert_$( basename $old )_contigs.err \
    "$CODEDIR/scripts/data_processing/convert_bigwig_contigs.py $old $new"
done < $TMPDIR/conservation_bw_urls.list
# GTEx eQTL
bsub \
  -q normal -J dl_GTEx_eQTL \
  -o $WRKDIR/LSF/logs/dl_GTEx_eQTL.log \
  -e $WRKDIR/LSF/logs/dl_GTEx_eQTL.err \
  "wget -P $TMPDIR/ https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz"
tar -xzf $TMPDIR/GTEx_Analysis_v7_eQTL.tar.gz -C $TMPDIR/
$CODEDIR/scripts/data_processing/gtex_eqtl_to_vcf.py \
  --gtex-tsv $TMPDIR/GTEx_Analysis_v7_eQTL/Pancreas.v7.signif_variant_gene_pairs.txt.gz \
  --gtex-tsv $TMPDIR/GTEx_Analysis_v7_eQTL/Colon_Sigmoid.v7.signif_variant_gene_pairs.txt.gz \
  --gtex-tsv $TMPDIR/GTEx_Analysis_v7_eQTL/Colon_Transverse.v7.signif_variant_gene_pairs.txt.gz \
  --gtex-tsv $TMPDIR/GTEx_Analysis_v7_eQTL/Lung.v7.signif_variant_gene_pairs.txt.gz \
  --header $WRKDIR/../refs/simple_hg19_header.somatic.vcf.gz \
| bcftools sort -O z -o $VEP_CACHE/GTEx_eQTLs.vcf.gz
tabix -p vcf -f $VEP_CACHE/GTEx_eQTLs.vcf.gz
# Gencode map of ENSG:symbol:ENST:CDS length
cat <<EOF > $WRKDIR/LSF/scripts/build_transcript_table.sh
#!/usr/bin/env bash
$CODEDIR/scripts/data_processing/build_gencode_table.py \
  --gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --outfile $TMPDIR/gencode.v19.annotation.transcript_info.unsorted.tsv \
  --no-header

sort -Vk3,3 -k2,2V -k1,1V -k4,4n \
  $TMPDIR/gencode.v19.annotation.transcript_info.unsorted.tsv \
| gzip -c > \
$WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz
EOF
chmod a+x $WRKDIR/LSF/scripts/build_transcript_table.sh
bsub \
  -q normal -J build_transcript_table \
  -o $WRKDIR/LSF/logs/build_transcript_table.log \
  -e $WRKDIR/LSF/logs/build_transcript_table.err \
  $WRKDIR/LSF/scripts/build_transcript_table.sh
# SpliceAI
for suf in gz gz.tbi; do
  ln -s \
    /data/talkowski/dg520/ref/spliceai/spliceai_scores.raw.snv.hg19.vcf.$suf \
    $VEP_CACHE/SpliceAI/spliceai_scores.raw.snv.hg37.vcf.$suf
done


######################
### VEP annotation ###
######################

### Write script for generic VEP function call with plugins and options
# export PERL5LIB=$VEP_PLUGINS/loftee:$PERL5LIB
cat <<EOF > $WRKDIR/LSF/scripts/run_VEP.sh
#!/usr/bin/env bash
set -eu -o pipefail

# Extract CNVs larger than 100bp to be spiked in after VEP
# (VEP will exclude these given --max_sv_size 100)
cnv_records="$TMPDIR/\$( basename \$1 | sed 's/vcf.gz/cnvs.vcf_records/g' )"
bcftools view \$1 --no-header \
  --include 'INFO/SVTYPE ~ "DEL\|DUP\|AMP"' \
  -O v -o \$cnv_records

# Annotate with VEP
$WRKDIR/../code/ensembl-vep/vep \
  --input_file \$1 \
  --format vcf \
  --output_file STDOUT \
  --vcf \
  --force_overwrite \
  --species homo_sapiens \
  --assembly GRCh37 \
  --max_sv_size 100 \
  --offline \
  --no_stats \
  --cache \
  --dir_cache $VEP_CACHE/ \
  --cache_version 108 \
  --buffer_size 2500 \
  --dir_plugins $VEP_PLUGINS/ \
  --fasta $TCGADIR/refs/GRCh37.fa \
  --minimal \
  --sift b \
  --polyphen b \
  --nearest gene \
  --distance 10000 \
  --numbers \
  --hgvs \
  --no_escape \
  --symbol \
  --canonical \
  --domains \
  --plugin UTRAnnotator,file=$VEP_CACHE/UTRAnnotator/uORF_5UTR_GRCh37_PUBLIC.txt \
  --plugin dbNSFP,$VEP_CACHE/dbNSFP4.3a_grch37.gz,FATHMM_score,MPC_score \
  --plugin CADD,$VEP_CACHE/CADD/whole_genome_SNVs.tsv.gz,$VEP_CACHE/CADD/InDels.tsv.gz \
  --plugin LoF,loftee_path:$VEP_PLUGINS/loftee,human_ancestor_fa:$VEP_CACHE/loftee/human_ancestor.fa.gz,conservation_file:$VEP_CACHE/loftee/phylocsf_gerp.sql \
  --plugin SpliceAI,snv=$VEP_CACHE/SpliceAI/spliceai_scores.raw.snv.hg37.vcf.gz,indel=$VEP_CACHE/SpliceAI/spliceai_scores.raw.indel.hg37.vcf.gz \
  --custom $VEP_CACHE/All_GRCh37_RS.bw,GERP,bigwig \
  --custom $VEP_CACHE/GRCh37.100way.phastCons.bw,phastCons,bigwig \
  --custom $VEP_CACHE/hg19.100way.phyloP100way.bw,phyloP,bigwig \
  --custom $VEP_CACHE/gnomad/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_POPMAX \
  --custom $VEP_CACHE/gnomad/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz,gnomADe,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_POPMAX \
  --custom $VEP_CACHE/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG \
  --custom $VEP_CACHE/cosmic/cmc.vcf.gz,COSMIC,vcf,exact,0,COSMIC_GENE,COSMIC_AA_CHANGE,COSMIC_GENE_TIER,COSMIC_MUT_SIG,COSMIC_MUT_FREQ \
  --custom $VEP_CACHE/gencode_promoters.bed.gz,promoters,bed,overlap,0 \
  --custom $VEP_CACHE/ABC_enhancers.bed.gz,ABC_enhancers,bed,overlap,0 \
  --custom $VEP_CACHE/CCR.bed.gz,CCR,bed,overlap,0 \
  --custom $VEP_CACHE/samocha_regional_constraint.bed.gz,ExAC_regional_constraint,bed,overlap,0 \
  --custom $VEP_CACHE/GTEx_eQTLs.vcf.gz,GTEx,vcf,exact,0,GTEx_eGene,GTEx_eQTL_beta,GTEx_eQTL_tissue \
| cat - \$cnv_records \
| bcftools sort -O z -o \$2
rm \$cnv_records
tabix -p vcf -f \$2
EOF
chmod a+x $WRKDIR/LSF/scripts/run_VEP.sh


### Annotate VCFs
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
  esac
  for subset in somatic_variants RAS_loci; do
    for suf in err log; do
      if [ -e $WRKDIR/LSF/logs/VEP_${cohort}_$subset.$suf ]; then
        rm $WRKDIR/LSF/logs/VEP_${cohort}_$subset.$suf
      fi
    done
    bsub -q big-multi -sla miket_sc -R "rusage[mem=16000]" -n 4 \
      -J VEP_${cohort}_$subset \
      -o $WRKDIR/LSF/logs/VEP_${cohort}_$subset.log \
      -e $WRKDIR/LSF/logs/VEP_${cohort}_$subset.err \
      "$WRKDIR/LSF/scripts/run_VEP.sh \
         $COHORTDIR/data/$cohort.$subset.vcf.gz \
         $COHORTDIR/data/$cohort.$subset.anno.vcf.gz"
  done
done


#########################
### Clean VEP outputs ###
#########################

### Write script for cleaning generic VEP output (generated above)
cat <<EOF > $WRKDIR/LSF/scripts/clean_VEP.sh
#!/usr/bin/env bash
set -eu -o pipefail

$CODEDIR/scripts/data_processing/cleanup_vep.py \
  --gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --transcript-info $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
  \$1 stdout \
| grep -ve "^##bcftools" | grep -ve "^##CADD_" | grep -ve "^##UTRAnnotator_" \
| grep -ve "^##LoF_" | grep -ve "^##SpliceAI_" | grep -ve "^##VEP-command-line" \
| bgzip -c \
> \$2
tabix -f -p vcf \$2
EOF
chmod a+x $WRKDIR/LSF/scripts/clean_VEP.sh


### Clean up VCFs
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      ;;
  esac
  for subset in somatic_variants RAS_loci; do
    for suf in err log; do
      if [ -e $WRKDIR/LSF/logs/VEP_cleanup_${cohort}_$subset.$suf ]; then
        rm $WRKDIR/LSF/logs/VEP_cleanup_${cohort}_$subset.$suf
      fi
    done
    bsub -q big -R "rusage[mem=16000]" -sla miket_sc \
      -J VEP_cleanup_${cohort}_$subset \
      -o $WRKDIR/LSF/logs/VEP_cleanup_${cohort}_$subset.log \
      -e $WRKDIR/LSF/logs/VEP_cleanup_${cohort}_$subset.err \
      "$WRKDIR/LSF/scripts/clean_VEP.sh \
         $COHORTDIR/data/$cohort.$subset.anno.vcf.gz \
         $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz"
  done
done


##############################
### Annotate in-sample AFs ###
##############################

### Write script to annotate in-cohort allele frequencies for all variants by population and cancer type
cat <<EOF > $WRKDIR/LSF/scripts/annotate_AFs.sh
#!/usr/bin/env bash
set -eu -o pipefail

$CODEDIR/scripts/data_processing/annotate_AFs.py \
  --sample-metadata \$2 \
  --sample-id-column \$3 \
  \$1 \
  \$4
tabix -f -p vcf \$4
EOF
chmod a+x $WRKDIR/LSF/scripts/annotate_AFs.sh


### Annotate AFs for all VCFs
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      sample_field=DONOR_ID
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      sample_field=PBP
      ;;
  esac
  for subset in somatic_variants RAS_loci; do
    for suf in err log; do
      if [ -e $WRKDIR/LSF/logs/annotate_AFs_${cohort}_$subset.$suf ]; then
        rm $WRKDIR/LSF/logs/annotate_AFs_${cohort}_$subset.$suf
      fi
    done
    bsub -q long -sla miket_sc \
      -J annotate_AFs_${cohort}_$subset \
      -o $WRKDIR/LSF/logs/annotate_AFs_${cohort}_$subset.log \
      -e $WRKDIR/LSF/logs/annotate_AFs_${cohort}_$subset.err \
      "$WRKDIR/LSF/scripts/annotate_AFs.sh \
         $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz \
         $COHORTDIR/data/sample_info/$cohort.ALL.sample_metadata.tsv.gz \
         $sample_field \
         $COHORTDIR/data/$cohort.$subset.anno.clean.wAF.vcf.gz"
  done
done


######################################
### Build simple genotype matrixes ###
######################################
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      sample_field=DONOR_ID
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      sample_field=PBP
      ;;
  esac
  for subset in somatic_variants RAS_loci; do
    for suf in err log; do
      if [ -e $WRKDIR/LSF/logs/make_dosage_${cohort}_$subset.$suf ]; then
        rm $WRKDIR/LSF/logs/make_dosage_${cohort}_$subset.$suf
      fi
    done
    bsub -q normal -sla miket_sc \
      -J make_dosage_${cohort}_$subset \
      -o $WRKDIR/LSF/logs/make_dosage_${cohort}_$subset.log \
      -e $WRKDIR/LSF/logs/make_dosage_${cohort}_$subset.err \
      "$CODEDIR/scripts/data_processing/vcf2dosage.py \
         $COHORTDIR/data/$cohort.$subset.anno.clean.vcf.gz - \
       | gzip -c > $COHORTDIR/data/$cohort.$subset.dosage.tsv.gz"
  done
done


################################################################
### Define sets of samples lacking DFCI/BWH Tier 1 mutations ###
################################################################
# Extract lists of samples
for cohort in TCGA PROFILE; do
  case $cohort in
    TCGA)
      COHORTDIR=$TCGADIR
      sample_field="donors"
      ;;
    PROFILE)
      COHORTDIR=$PROFILEDIR
      sample_field="samples"
      ;;
  esac
  $CODEDIR/scripts/data_processing/define_control_samples.py \
    --vcf $COHORTDIR/data/$cohort.somatic_variants.anno.clean.vcf.gz \
    --criteria $CODEDIR/refs/RAS_control_sample_criteria.json \
    --regions $CODEDIR/refs/RAS_loci.GRCh37.bed.gz \
    --eligible-samples $COHORTDIR/data/sample_info/$cohort.ALL.$sample_field.list \
    --exclude-samples $COHORTDIR/data/sample_info/$cohort.ALL.$sample_field.missing_somatic.list \
    --outfile $COHORTDIR/data/sample_info/$cohort.ALL.eligible_controls.list
done
# Summarize as table
for cancer in PDAC CRAD LUAD SKCM; do
  echo $cancer
  for cohort in TCGA PROFILE; do
    case $cohort in
      TCGA)
        COHORTDIR=$TCGADIR
        sample_field="donors"
        ;;
      PROFILE)
        COHORTDIR=$PROFILEDIR
        sample_field="samples"
        ;;
    esac
    elig_samps=$COHORTDIR/data/sample_info/$cohort.$cancer.$sample_field.list
    fgrep -wf $elig_samps \
      $COHORTDIR/data/sample_info/$cohort.ALL.eligible_controls.list \
    | wc -l | addcom
  done | paste -s -
done | paste - -
