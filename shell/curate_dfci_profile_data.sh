#!/usr/bin/env bash

################################
#    EGFR Modifiers Project    #
################################

# Copyright (c) 2023-Present Ryan L. Collins, and the Gusev/Van Allen Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Curate data from DFCI-PROFILE cohort for EGFR modifiers study

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export BASEDIR=/data/gusev/PROFILE
export GTDIR=$BASEDIR/2020_2022_combined/IMPUTE_HQ
export CLINDIR=$BASEDIR/CLINICAL
export WRKDIR=/data/gusev/USERS/rlc47/PROFILE
export CODEDIR=$WRKDIR/../EGFR_modifier_analysis/code/ras_modifiers
cd $WRKDIR


### NOTE: EXPECTS FULL RAS MODIFIER DATA CURATION TO BE ALREADY COMPLETE 
### SEE VERSION OF THIS FILE IN THE MAIN BRANCH


### Subset main patient metadata .tsv to just LUAD samples for convenience
all_sample_meta=$WRKDIR/data/sample_info/PROFILE.ALL.sample_metadata.tsv.gz
zcat $all_sample_meta \
| fgrep -wf $WRKDIR/data/sample_info/PROFILE.LUAD.samples.list \
| cat <( zcat $all_sample_meta | head -n1 ) - \
| gzip -c \
> $( echo $all_sample_meta | sed 's/\.ALL\./.LUAD./g' )


# Extract samples & loci of interest
while read contig start end gene; do
  cat << EOF > $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools view \
  -O z -o $WRKDIR/data/PROFILE.$gene.vcf.gz \
  --min-ac 1 \
  --samples-file $WRKDIR/data/sample_info/PROFILE.LUAD.samples.list \
  --regions "$contig:${start}-$end" \
  $GTDIR/PROFILE_COMB.$contig.HQ.vcf.gz
tabix -p vcf -f $WRKDIR/data/PROFILE.$gene.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
  rm $WRKDIR/LSF/logs/extract_${gene}_variants.*
  bsub \
    -q normal -R 'rusage[mem=6000]' -n 2 -J PROFILE_extract_${gene}_variants \
    -o $WRKDIR/LSF/logs/extract_${gene}_variants.log \
    -e $WRKDIR/LSF/logs/extract_${gene}_variants.err \
    $WRKDIR/LSF/scripts/extract_${gene}_variants.sh
done < <( zcat $CODEDIR/refs/EGFR_loci.GRCh37.bed.gz | fgrep -v "#" )
# Merge VCFs for each gene into a single VCF and index the merged VCF
bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
  --header-lines <( zcat $WRKDIR/refs/simple_hg19_header.vcf.gz | fgrep -w SVTYPE ) \
  -O z -o $WRKDIR/data/PROFILE.EGFR_loci.vcf.gz \
  $WRKDIR/data/PROFILE.EGFR.vcf.gz
tabix -p vcf -f $WRKDIR/data/PROFILE.EGFR_loci.vcf.gz


# Summarize somatic variant status by gene & cancer type
n_samp=$( fgrep -wvf \
            $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic.list \
            ${WRKDIR}/data/sample_info/PROFILE.LUAD.samples.list | wc -l )
while read gene; do
  zcat $WRKDIR/data/PROFILE.somatic_variants.tsv.gz \
  | fgrep -wf ${WRKDIR}/data/sample_info/PROFILE.LUAD.samples.list \
  | fgrep -wvf $WRKDIR/data/sample_info/PROFILE.ALL.samples.missing_somatic.list \
  | awk -v gene=$gene -v FS="\t" '{ if ($9==gene && ($15~/AMP|DEL|Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|Splice_Site/)) print $7 }' \
  | sort | uniq | wc -l | awk -v n=$n_samp '{ print $1/n }'
done < <( zcat $CODEDIR/refs/EGFR_loci.GRCh37.bed.gz | fgrep -v "#" | cut -f4 )

