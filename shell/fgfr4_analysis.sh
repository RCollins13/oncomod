#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Focused investigation of germline FGFR4 association with KRAS tier 1 mutations

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export TCGADIR=/data/gusev/USERS/rlc47/TCGA
export PROFILEDIR=/data/gusev/USERS/rlc47/PROFILE
export HMFDIR=/data/gusev/USERS/rlc47/HMF
export WRKDIR=/data/gusev/USERS/rlc47/RAS_modifier_analysis
export CODEDIR=$WRKDIR/../code/oncomod
cd $WRKDIR


### Prep directory structure
for dir in $WRKDIR/data/FGFR4 $WRKDIR/results/FGFR4; do
  if ! [ -e $dir ]; then mkdir $WRKDIR/data/FGFR4; fi
done


### Compute single-variant association tests for all FGFR4 coding variants vs. all KRAS somatic endpoints
# Make input germline test sets
for cohort in TCGA PROFILE HMF; do
  zcat $WRKDIR/data/variant_sets/$cohort.germline.collapsed_coding_csqs.tsv.gz \
  | awk -v FS="\t" -v OFS="\t" '{ if ($2=="FGFR4") print $1, $NF }' \
  | sort -Vk1,1
done \
| $CODEDIR/scripts/germline_somatic_assoc/simple_collapse_variant_sets.py \
  --tsv-out $WRKDIR/data/FGFR4/FGFR4.coding_variants.sets.tsv
# Run pooled association tests for each FGFR4 coding variant
cat << EOF > $WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.sh
#!/usr/bin/env bash
$CODEDIR/scripts/germline_somatic_assoc/germline_somatic_assoc.pooled.R \
--sample-metadata $TCGADIR/data/sample_info/TCGA.ALL.sample_metadata.tsv.gz \
--somatic-ad $TCGADIR/data/TCGA.somatic_variants.dosage.tsv.gz \
--germline-ad $TCGADIR/data/TCGA.RAS_loci.dosage.tsv.gz \
--name TCGA \
--germline-gq $TCGADIR/data/TCGA.RAS_loci.GQ.tsv.gz \
--eligible-controls $TCGADIR/data/sample_info/TCGA.ALL.eligible_controls.list \
--sample-metadata $PROFILEDIR/data/sample_info/PROFILE.ALL.sample_metadata.tsv.gz \
--somatic-ad $PROFILEDIR/data/PROFILE.somatic_variants.dosage.tsv.gz \
--germline-ad $PROFILEDIR/data/PROFILE.RAS_loci.dosage.tsv.gz \
--name DFCI \
--germline-gq $PROFILEDIR/data/PROFILE.RAS_loci.GQ.tsv.gz \
--eligible-controls $PROFILEDIR/data/sample_info/PROFILE.ALL.eligible_controls.list \
--sample-metadata $HMFDIR/data/sample_info/HMF.ALL.sample_metadata.tsv.gz \
--somatic-ad $HMFDIR/data/HMF.somatic_variants.dosage.tsv.gz \
--germline-ad $HMFDIR/data/HMF.RAS_loci.dosage.tsv.gz \
--name HMF \
--germline-gq $HMFDIR/data/HMF.RAS_loci.GQ.tsv.gz \
--eligible-controls $HMFDIR/data/sample_info/HMF.ALL.eligible_controls.list \
--cancer-type CRAD \
--somatic-variant-sets $WRKDIR/data/variant_sets/test_sets/CRAD.KRAS.somatic_endpoints.tsv \
--germline-variant-sets $WRKDIR/data/FGFR4/FGFR4.coding_variants.sets.tsv \
--outfile $WRKDIR/results/FGFR4/pooled.CRAD.FGFR4_allelic_series.sumstats.tsv
gzip -f $WRKDIR/results/FGFR4/pooled.CRAD.FGFR4_allelic_series.sumstats.tsv
EOF
chmod a+x $WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.sh
for suf in err log; do
  logfile=$WRKDIR/LSF/logs/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.$suf
  if [ -e $logfile ]; then rm $logfile; fi
done
bsub -q big-multi -sla miket_sc -R "rusage[mem=24000]" -n 4 \
  -J germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series \
  -o $WRKDIR/LSF/logs/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.log \
  -e $WRKDIR/LSF/logs/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.err \
  "$WRKDIR/LSF/scripts/germline_somatic_assoc_pooled_CRAD.FGFR4_allelic_series.sh"


