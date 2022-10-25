#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Locus refinement for targeted KRAS/NRAS/HRAS analyses

# Note: intended to be executed on the Broad cluster


### Set local parameters
export BASEDIR=/broad/VanAllenLab/xchip/cga_home/rcollins
export WRKDIR=/broad/VanAllenLab/xchip/cga_home/rcollins/ras_modifiers
export CODEDIR=$WRKDIR/code/ras_modifiers
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in refs data code; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Write simple reference files
echo -e "NRAS\nHRAS\nKRAS" > $WRKDIR/refs/genes.list


### Extract gene bodies from MANE-Select GTF
wget -O - https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz \
| gunzip -c | sort -Vk1,1 -k4,4n -k5,5n | bgzip -c \
> $WRKDIR/refs/MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz
tabix -p gff -f $WRKDIR/refs/MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz
export GTF=$WRKDIR/refs/MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz
zcat $GTF | fgrep -wf $WRKDIR/refs/genes.list | fgrep -w MANE_Select \
| awk -v OFS="\t" '{ if ($3=="transcript") print $1, $4, $5, $16 }' \
| tr -d '";' | bgzip -c > $WRKDIR/refs/genes.bed.gz


### Extract ENSG-to-symbol mappings for genes of interest
zcat $GTF | fgrep -wf $WRKDIR/refs/genes.list \
| awk -v OFS="\t" '{ if ($3=="gene") print $10, $14 }' \
| tr -d '";' > $WRKDIR/refs/ENSG_to_symbol.tsv


### Find most distant eQTLs from GTEx v8 (any tissue)
wget -O - https://storage.googleapis.com/gtex_analysis_v8/multi_tissue_qtl_data/GTEx_Analysis_v8.metasoft.txt.gz \
| gunzip -c | head -n1 > $WRKDIR/data/GTEx_Analysis_v8.metasoft.RAS_genes.tsv
wget -O - https://storage.googleapis.com/gtex_analysis_v8/multi_tissue_qtl_data/GTEx_Analysis_v8.metasoft.txt.gz \
| gunzip -c | fgrep -f <( cut -f1 -d \. $WRKDIR/refs/ENSG_to_symbol.tsv ) \
>> $WRKDIR/data/GTEx_Analysis_v8.metasoft.RAS_genes.tsv
gzip -f $WRKDIR/data/GTEx_Analysis_v8.metasoft.RAS_genes.tsv
$CODEDIR/scripts/data_processing/preprocess_eQTLs.py \
  --ensg-map $WRKDIR/refs/ENSG_to_symbol.tsv \
  --outfile $WRKDIR/data/RAS_eQTLs.tsv.gz \
  --verbose-outfile $WRKDIR/data/RAS_eQTLs.verbose.tsv.gz \
  $WRKDIR/data/GTEx_Analysis_v8.metasoft.RAS_genes.tsv.gz
while read gene; do
  zcat $WRKDIR/data/RAS_eQTLs.tsv.gz \
  | fgrep -w $gene | cut -f2 | sort -nrk1,1 | sed -e 1b -e '$!d' | paste -s -
done < $WRKDIR/refs/genes.list