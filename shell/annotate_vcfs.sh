#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Annotate all VCFs prior to analysis

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export CODEDIR=$WRKDIR/../code/ras_modifiers
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in LSF LSF/logs LSF/scripts; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Download & install VEP and associated plugins
# Prep directory structure
if ! [ -e $WRKDIR/../refs/vep_cache ]; then
  mkdir $WRKDIR/../refs/vep_cache
fi
if ! [ -e $WRKDIR/../code/vep_plugins ]; then
  mkdir $WRKDIR/../code/vep_plugins
fi
cd $WRKDIR/../code/ && \
git clone https://github.com/Ensembl/ensembl-vep.git && \
cd ensembl-vep && \
perl INSTALL.pl \
  --NO_HTSLIB \
  --CONVERT \
  --CACHEDIR $WRKDIR/../refs/vep_cache/ \
  --PLUGINSDIR $WRKDIR/../code/vep_plugins/ && \
cd $WRKDIR
# Download GRCh37 & GRCh38 indexed VEP caches
for ref in GRCh37 GRCh38; do
  bsub \
    -q normal -J ${ref}_VEP_cache_download \
    -o $WRKDIR/LSF/logs/${ref}_VEP_cache_download.log \
    -e $WRKDIR/LSF/logs/${ref}_VEP_cache_download.err \
    "cd $WRKDIR/../refs/vep_cache && \
     wget https://ftp.ensembl.org/pub/release-108/variation/indexed_vep_cache/homo_sapiens_vep_108_${ref}.tar.gz && \
     tar -xzf homo_sapiens_vep_108_${ref}.tar.gz"
done
# Test VEP installation
cd $WRKDIR/../code/ensembl-vep && \
./vep \
  -i examples/homo_sapiens_GRCh37.vcf \
  --cache \
  --dir_cache $WRKDIR/../refs/vep_cache/ \
  --cache_version 108 \
  --offline && \
cd $WRKDIR


### Curate custom datasets for use by VEP
# Prep subdirectories
for subdir in UTRAnnotator gnomad loftee CADD clinvar cosmic; do
  if ! [ -e $WRKDIR/../refs/vep_cache/$subdir ]; then
    mkdir $WRKDIR/../refs/vep_cache/$subdir
  fi
done
# dbNSFP
bsub \
  -q normal -J dbNSFP_cache_download \
  -o $WRKDIR/LSF/logs/dbNSFP_cache_download.log \
  -e $WRKDIR/LSF/logs/dbNSFP_cache_download.err \
  "cd $WRKDIR/../refs/vep_cache && \
   wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.3a.zip && \
   unzip dbNSFP4.3a.zip"
zcat $WRKDIR/../refs/vep_cache/dbNSFP4.3a_variant.chr1.gz \
| head -n1 > $TMPDIR/dbNSFP_header
cat <<EOF > $WRKDIR/LSF/scripts/prep_dbNSFP.grch38.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
zgrep -h -v ^#chr $WRKDIR/../refs/vep_cache/dbNSFP4.3a_variant.chr* \
| sort -T $TMPDIR -k1,1 -k2,2n - | cat $TMPDIR/dbNSFP_header - \
| bgzip -c > $WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch38.gz
tabix -s 1 -b 2 -e 2 -f $WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch38.gz
EOF
cat <<EOF > $WRKDIR/LSF/scripts/prep_dbNSFP.grch37.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
zgrep -h -v ^#chr $WRKDIR/../refs/vep_cache/dbNSFP4.3a_variant.chr* \
| awk '\$8 != "." ' | sort -T $TMPDIR -k8,8 -k9,9n - | cat $TMPDIR/dbNSFP_header - \
| bgzip -c > $WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch37.gz
tabix -s 8 -b 9 -e 9 -f $WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch37.gz
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
for file in $WRKDIR/../code/vep_plugins/loftee/*; do
  ln -fs $( echo -e $file ) $WRKDIR/../code/vep_plugins/$( basename $file )
done
for suffix in gz gz.fai gz.gzi; do
  bsub \
    -q normal -J dl_loftee_ancestor_fa$suffix \
    -o $WRKDIR/LSF/logs/dl_loftee_ancestor_fa$suffix.log \
    -e $WRKDIR/LSF/logs/dl_loftee_ancestor_fa$suffix.err \
    "cd $WRKDIR/../refs/vep_cache/loftee && \
     wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.$suffix"
done
bsub \
  -q normal -J dl_loftee_conservation_db \
  -o $WRKDIR/LSF/logs/dl_loftee_conservation_db.log \
  -e $WRKDIR/LSF/logs/dl_loftee_conservation_db.err \
  "cd $WRKDIR/../refs/vep_cache/loftee && \
   wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz && \
   gunzip phylocsf_gerp.sql.gz"
# CADD
for prefix in whole_genome_SNVs InDels; do
  for suffix in gz gz.tbi; do
    bsub \
      -q normal -J dl_cadd_$prefix.tsv.$suffix \
      -o $WRKDIR/LSF/logs/dl_cadd_$prefix.tsv.$suffix.log \
      -e $WRKDIR/LSF/logs/dl_cadd_$prefix.tsv.$suffix.err \
      "cd $WRKDIR/../refs/vep_cache/CADD && \
       wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/$prefix.tsv.$suffix"
  done
done
# ClinVar
wget \
  -P $WRKDIR/../refs/vep_cache/clinvar/ \
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
   | bcftools sort -O z -o $WRKDIR/../refs/vep_cache/cosmic/cmc.vcf.gz && \
   tabix -p vcf -f $WRKDIR/../refs/vep_cache/cosmic/cmc.vcf.gz"
# ABC
wget -P $TMPDIR/ \
  ftp://ftp.broadinstitute.org/outgoing/lincRNA/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz
zcat $TMPDIR/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz
zcat $TMPDIR/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz \
| awk '{ print $NF }' | sort | uniq -c
# Tissues of interest: pancreas-Roadmap Panc1-ENCODE body_of_pancreas-ENCODE
# sigmoid_colon-ENCODE transverse_colon-ENCODE keratinocyte-Roadmap fibroblast_of_dermis-Roadmap fibroblast_of_arm-ENCODE
# ftp://ftp.broadinstitute.org/outgoing/lincRNA/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz
# gnomAD
for subset in exomes genomes; do
  bsub \
    -q normal -J gnomad_${subset}_cache_download \
    -o $WRKDIR/LSF/logs/gnomad_${subset}_cache_download.log \
    -e $WRKDIR/LSF/logs/gnomad_${subset}_cache_download.err \
    "cd $WRKDIR/../refs/vep_cache/gnomad && \
     wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.$subset.r2.0.1.sites.noVEP.vcf.gz && \
     wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.$subset.r2.0.1.sites.noVEP.vcf.gz.tbi"
done
# Gene promoters (defined as +2kb upstream of TSS)
$CODEDIR/scripts/data_processing/make_promoters.py \
  --gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --coding-only \
| sort -Vk1,1 -k2,2n -k3,3n | bgzip -c \
> $WRKDIR/../refs/vep_cache/gencode_promoters.bed.gz
tabix -p bed -f $WRKDIR/../refs/vep_cache/gencode_promoters.bed.gz


### Build script for generic VEP function call with plugins and options
# export PERL5LIB=$WRKDIR/../code/vep_plugins/loftee:$PERL5LIB
annotate_vcf () {
  $WRKDIR/../code/ensembl-vep/vep \
    --input_file $1 \
    --format vcf \
    --output_file $2 \
    --vcf \
    --compress_output bgzip \
    --force_overwrite \
    --species homo_sapiens \
    --assembly GRCh37 \
    --offline \
    --no_stats \
    --fork 2 \
    --cache \
    --dir_cache $WRKDIR/../refs/vep_cache/ \
    --cache_version 108 \
    --dir_plugins $WRKDIR/../code/vep_plugins/ \
    --fasta $TCGADIR/refs/GRCh37.fa \
    --minimal \
    --sift b \
    --polyphen b \
    --nearest gene \
    --regulatory \
    --cell_type pancreas,endocrine_pancreas,sigmoid_colon,foreskin_melanocyte_1,foreskin_melanocyte_2 \
    --numbers \
    --hgvs \
    --symbol \
    --canonical \
    --domains \
    --plugin UTRAnnotator,file=$WRKDIR/../refs/vep_cache/UTRAnnotator/uORF_5UTR_GRCh37_PUBLIC.txt \
    --plugin dbNSFP,$WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch37.gz,SIFT_score,Polyphen2_HDIV_score,FATHMM_score,MPC_score,GERP++_RS \
    --plugin CADD,$WRKDIR/../refs/vep_cache/CADD/whole_genome_SNVs.tsv.gz,$WRKDIR/../refs/vep_cache/CADD/InDels.tsv.gz \
    --plugin LoF,loftee_path:$WRKDIR/../code/vep_plugins/loftee,human_ancestor_fa:$WRKDIR/../refs/vep_cache/loftee/human_ancestor.fa.gz,conservation_file:$WRKDIR/../refs/vep_cache/loftee/phylocsf_gerp.sql \
    --custom $WRKDIR/../refs/vep_cache/gnomad/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_POPMAX \
    --custom $WRKDIR/../refs/vep_cache/gnomad/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz,gnomADe,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_POPMAX \
    --custom $WRKDIR/../refs/vep_cache/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --custom $WRKDIR/../refs/vep_cache/cosmic/cmc.vcf.gz,COSMIC,vcf,exact,0,COSMIC_GENE,COSMIC_AA_CHANGE,COSMIC_GENE_TIER,COSMIC_MUT_SIG,COSMIC_MUT_FREQ \
    --custom $WRKDIR/../refs/vep_cache/gencode_promoters.bed.gz,promoters,bed,exact,0
}
    
# DEV:
annotate_vcf $TMPDIR/test.vcf.gz $TMPDIR/test.anno.vcf.gz


### Annotate TCGA VCFs
annotate_vcf \
  $TCGADIR/data/TCGA.RAS_loci.vcf.gz \
  $TCGADIR/data/TCGA.RAS_loci.anno.vcf.gz

### Annotate PROFILE VCFs
# TODO: implement this

