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
for SUBDIR in LSF LSF/logs; do
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
zgrep -h -v ^#chr $WRKDIR/../refs/vep_cache/dbNSFP4.3a_variant.chr* \
| sort -T $TMPDIR -k1,1 -k2,2n - | cat $TMPDIR/dbNSFP_header - \
| bgzip -c > $WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch38.gz
tabix -s 1 -b 2 -e 2 $WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch38.gz
zgrep -h -v ^#chr $WRKDIR/../refs/vep_cache/dbNSFP4.3a_variant.chr* \
| awk '$8 != "." ' | sort -T $TMPDIR -k8,8 -k9,9n - | cat $TMPDIR/dbNSFP_header - \
| bgzip -c > $WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch37.gz
tabix -s 8 -b 9 -e 9 $WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch37.gz
# LOFTEE
# TODO: implement this
# ABC
# TODO: implement this
# ftp://ftp.broadinstitute.org/outgoing/lincRNA/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz
# gnomAD
if ! [ -e $WRKDIR/../refs/vep_cache/gnomad ]; then
  mkdir $WRKDIR/../refs/vep_cache/gnomad
fi
for subset in exomes genomes; do
  bsub \
    -q normal -J gnomad_${subset}_cache_download \
    -o $WRKDIR/LSF/logs/gnomad_${subset}_cache_download.log \
    -e $WRKDIR/LSF/logs/gnomad_${subset}_cache_download.err \
    "cd $WRKDIR/../refs/vep_cache/gnomad && \
     wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.$subset.r2.0.1.sites.noVEP.vcf.gz && \
     wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.$subset.r2.0.1.sites.noVEP.vcf.gz.tbi"
done


### Build generic VEP function call with plugins and options
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
    --af_gnomadg \
    --plugin UTRAnnotator,file=$WRKDIR/../refs/vep_cache/UTRAnnotator/uORF_5UTR_GRCh37_PUBLIC.txt
}
# --plugin dbNSFP,$WRKDIR/../refs/vep_cache/dbNSFP4.3a_grch37.gz,LRT_score,GERP++_RS


### Annotate TCGA VCFs
annotate_vcf \
  $TCGADIR/data/TCGA.RAS_loci.vcf.gz \
  $TCGADIR/data/TCGA.RAS_loci.anno.vcf.gz

### Annotate PROFILE VCFs
# TODO: implement this

