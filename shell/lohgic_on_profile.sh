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
for dir in $WRKDIR/LOHGIC $WRKDIR/LOHGIC/AllFIT $WRKDIR/LOHGIC/AllFIT/inputs \
           $WRKDIR/LOHGIC/AllFIT/outputs; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done


### Install All-FIT for purity estimation
cd $CODEDIR && \
git clone https://github.com/KhiabanianLab/All-FIT.git && \
cd All-FIT && \
chmod a+x All-FIT.py && \
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

# Step 2.1: run All-FIT from the outputs of step 1
cat << EOF > $WRKDIR/LSF/scripts/AllFIT.sh
#!/usr/bin/env bash

. /PHShome/rlc47/.bashrc
cd $WRKDIR

ID=\$1
OUTDIR=$WRKDIR/LOHGIC/AllFIT/outputs/\$ID

if [ -e \$OUTDIR ]; then
  rm -rf \$OUTDIR
fi
$CODEDIR/All-FIT/All-FIT.py \
  -i $WRKDIR/LOHGIC/AllFIT/inputs/\$ID.AllFIT_input.tsv \
  -d \$OUTDIR \
  -o \$ID \
  -t all
EOF
chmod a+x $WRKDIR/LSF/scripts/AllFIT.sh
while read infile; do
  ID=$( basename $infile | sed 's/\.AllFIT_input\.tsv//g' )
  bsub \
    -q vshort -sla miket_sc -J AllFIT_$ID \
    -o $WRKDIR/LSF/logs/AllFIT_$ID.log \
    -e $WRKDIR/LSF/logs/AllFIT_$ID.err \
    "$WRKDIR/LSF/scripts/AllFIT.sh $ID"
done < <( find $WRKDIR/LOHGIC/AllFIT/inputs/ -name "*.AllFIT_input.tsv" )

# Step 2.2: get list of samples included in AllFIT run and compile purity from 
# OncoPanel report for comparison
find $WRKDIR/LOHGIC/AllFIT/inputs/ \
  -name "*.AllFIT_input.tsv" \
| xargs -I {} basename {} \
| sed 's/\.AllFIT_input\.tsv//g' \
> $WRKDIR/LOHGIC/AllFIT/samples.list
cat << EOF > $TMPDIR/parse_oncdrs.R
#!/usr/bin/env Rscript
options(stringsAsFactors=F)
x <- read.table("$BASEDIR/CLINICAL/OncDRS/ALL_2022_11/GENOMIC_SPECIMEN.csv",
                header=T, sep=",")[, c(5, 2, 14, 12)]
s <- read.table("$WRKDIR/LOHGIC/AllFIT/samples.list", header=F)[, 1]
write.table(x[which(x\$SAMPLE_ACCESSION_NBR %in% s), ],
            "$WRKDIR/LOHGIC/AllFIT/sample_metadata.tsv",
            sep="\t", col.names=T, row.names=F, quote=F)
EOF
Rscript $TMPDIR/parse_oncdrs.R

# Step 3: collect All-FIT purity estimates for each sample
find $WRKDIR/LOHGIC/AllFIT/outputs/ -name "*.txt" \
| xargs -I {} tail -n1 {} \
| sed 's/,/\t/g' \
| awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, $NF }' \
| cat <( echo -e "SAMPLE_ACCESSION_NBR\tPURITY\tCI95_LOWER\tCI95_UPPER" ) - \
>  $WRKDIR/LOHGIC/AllFIT/PROFILE.AllFIT_purity_estimates.tsv

# Step 4: prepare data for LOHGIC
# LOHGIC requires .tsv as input with at least ploidy, VAF, total depth, and sample purity
# Optionally, VAF CI and purity CI can be provided as columns 5 & 6
cat << EOF > $WRKDIR/LSF/scripts/LOHGIC.sh
#!/usr/bin/env bash

. /PHShome/rlc47/.bashrc
cd $WRKDIR

module load matlab/default




