#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Infer germline coding variants from OncoPanel data for PROFILE cohort

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export BASEDIR=/data/gusev/PROFILE
export WRKDIR=/data/gusev/USERS/rlc47/PROFILE
export CODEDIR=$WRKDIR/../code
cd $WRKDIR


### Prepare directory tree
for dir in $WRKDIR/LOHGIC $WRKDIR/LOHGIC/AllFIT $WRKDIR/LOHGIC/AllFIT/inputs; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done


### Install All-FIT for purity estimation
cd $CODEDIR && \
git clone https://github.com/KhiabanianLab/All-FIT.git && \
cd All-FIT && \
ln -s `pwd`/All-FIT.py $CODEDIR/bin/All-FIT.py && \
cd $WRKDIR


##############################################################
### Run All-FIT to estimate purity for all PROFILE samples ###
##############################################################
# All-FIT must be called once per sample and expects a four-column .tsv as input
# Columns are: uniqueID, variant allele frequencies (*100), depth, and estimated ploidy
# See https://github.com/KhiabanianLab/All-FIT

# Step 1: unify mutation and copy number data for each biopsy to extract
# information required by All-FIT
$CODEDIR/ras_modifiers/scripts/lohgic/get_allfit_inputs.py \
  --mutations $BASEDIR/CLINICAL/OncDRS/ALL_2022_11/GENOMIC_MUTATION_RESULTS.csv \
  --cnas $BASEDIR/CLINICAL/OncDRS/ALL_2022_11/GENOMIC_CNV_RESULTS.csv \
  --outdir $WRKDIR/LOHGIC/AllFIT/inputs

# Step 2: run All-FIT from the outputs of step 1
# cat << EOF > $WRKDIR/LSF/scripts/AllFIT.sh
# #!/usr/bin/env bash
# . /PHShome/rlc47/.bashrc
# cd $WRKDIR

# EOF
# while read infile; do
#   sid=$( basename $infile | sed 's/\.AllFIT_input\.tsv//g' )
#   echo $sid
#   cat << EOF > $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
# #!/usr/bin/env bash
# . /PHShome/rlc47/.bashrc
# cd $WRKDIR
# bcftools view \
#   -O z -o $WRKDIR/data/PROFILE.$gene.vcf.gz \
#   --min-ac 1 \
#   --samples-file $WRKDIR/data/sample_info/PROFILE.ALL.samples.list \
#   --regions "$contig:${start}-$end" \
#   $GTDIR/PROFILE_COMB.$contig.HQ.vcf.gz
# tabix -p vcf -f $WRKDIR/data/PROFILE.$gene.vcf.gz
# EOF
#   chmod a+x $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
#   rm $WRKDIR/LSF/logs/extract_${gene}_variants.*
#   bsub \
#     -q normal -R 'rusage[mem=6000]' -n 2 -J PROFILE_extract_${gene}_variants \
#     -o $WRKDIR/LSF/logs/extract_${gene}_variants.log \
#     -e $WRKDIR/LSF/logs/extract_${gene}_variants.err \
#     $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
# done < <( zcat $CODEDIR/refs/RAS_loci.GRCh37.bed.gz | fgrep -v "#" )



