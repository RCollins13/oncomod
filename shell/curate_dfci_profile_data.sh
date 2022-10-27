#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2022-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Locus refinement for targeted KRAS/NRAS/HRAS analyses

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export BASEDIR=/data/gusev/PROFILE
export GTDIR=$BASEDIR/2020_2022_combined/IMPUTE_HQ
export CLINDIR=$BASEDIR/CLINICAL
export WRKDIR=/data/gusev/USERS/rlc47/PROFILE
export CODEDIR=$WRKDIR/code/ras_modifiers
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/sample_info; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Curate clinical information for patients of interest
# Get list of all IDs present in VCF
tabix -H $GTDIR/PROFILE_COMB.22.HQ.vcf.gz \
| fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
> $WRKDIR/data/sample_info/vcf.ids.list
# Get list of patients from cancer types of interest
# $CODEDIR/scripts/data_processing/preprocess_dfci_profile_ehr.py \
$TMPDIR/preprocess_dfci_profile_ehr.py \
  --id-map-tsv $CLINDIR/PROFILE_MRN_BL_PANEL.PBP.tab \
  --dx-csv $CLINDIR/OncDRS/ALL_2021_11/CANCER_DIAGNOSIS_CAREG.csv.gz \
  --ancestry-csv $CLINDIR/PROFILE_2022_ANCESTRY.csv.gz \
  --hx-csv $CLINDIR/OncDRS/ALL_2021_11/HEALTH_HISTORY.csv.gz \
  --survival-csv $CLINDIR/OncDRS/ALL_2021_11/PT_INFO_STATUS_REGISTRATION.csv.gz \
  --out-prefix $WRKDIR/data/sample_info/ \
  --vcf-ids $WRKDIR/data/sample_info/vcf.ids.list


### Subset VCFs to patients of interest and RAS loci
# Note: may need to lift over to hg38?


### Curate somatic data for patients of interest


# The imputed germline VCFs are here:
# /data/gusev/PROFILE/2020_2022_combined/IMPUTE_HQ
# And the samples use “PBP” identifiers
#  
# The identifier mapping is here:
# /data/gusev/PROFILE/CLINICAL/PROFILE_MRN_BL_PANEL.PBP.tab
# The PBP ids are for the imputed data, the BL ids are for the tumor specimen, and the MRN is the medical number for the patient. All but the PBP ids are sensitive data so make sure not to share them.
#  
# If you need it, genetic ancestry is here:
# /data/gusev/PROFILE/CLINICAL/PROFILE_2022_ANCESTRY.csv.gz 
#  
# All of the EHR data is here:
# /data/gusev/PROFILE/CLINICAL/OncDRS/ALL_2021_11
# (also sensitive data so please don’t share), the GENOMIC_SPECIMEN file contains information on the tumor biopsied, the *MUTATION* tables list the SNVs and the *CNV* tables list the CNVs, the CANCER_DIAGNOSIS_CAREG table contains cancer stage and other information from the cancer registry (on a subset of patients) which may also be useful to stratify by stage. The schema for these tables is described here: https://wiki.dfci.harvard.edu:8443/oncdata/latest
