#!/usr/bin/env bash

################################
#    EGFR Modifiers Project    #
################################

# Copyright (c) 2023-Present Ryan L. Collins, Jackie LoPiccolo, and the Gusev/Van Allen Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Locus refinement for targeted EGFR analyses

# Note: intended to be executed on the Broad cluster


### Set local parameters
export BASEDIR=/broad/VanAllenLab/xchip/cga_home/rcollins
export WRKDIR=/broad/VanAllenLab/xchip/cga_home/rcollins/egfr_modifiers
export CODEDIR=$WRKDIR/code/ras_modifiers
cd $CODEDIR && \
git checkout EGFR && \
git pull && \
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in refs data code UGER UGER/logs; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Write simple reference files
echo -e "EGFR" > $WRKDIR/refs/genes.list


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
| gunzip -c | head -n1 > $WRKDIR/data/GTEx_Analysis_v8.metasoft.EGFR.tsv
wget -O - https://storage.googleapis.com/gtex_analysis_v8/multi_tissue_qtl_data/GTEx_Analysis_v8.metasoft.txt.gz \
| gunzip -c | fgrep -f <( cut -f1 -d \. $WRKDIR/refs/ENSG_to_symbol.tsv ) \
>> $WRKDIR/data/GTEx_Analysis_v8.metasoft.EGFR.tsv
gzip -f $WRKDIR/data/GTEx_Analysis_v8.metasoft.EGFR.tsv
$CODEDIR/scripts/data_processing/preprocess_eQTLs.py \
  --ensg-map $WRKDIR/refs/ENSG_to_symbol.tsv \
  --outfile $WRKDIR/data/EGFR_eQTLs.tsv.gz \
  --verbose-outfile $WRKDIR/data/EGFR_eQTLs.verbose.tsv.gz \
  $WRKDIR/data/GTEx_Analysis_v8.metasoft.EGFR.tsv.gz
while read gene; do
  zcat $WRKDIR/data/EGFR_eQTLs.tsv.gz \
  | fgrep -w $gene | cut -f2 | sort -nk1,1 | sed -e 1b -e '$!d' | paste -s -
done < $WRKDIR/refs/genes.list


### Find most permissive TAD boundaries per gene in any tissue
# Note: need to have downloaded & unpacked data from Schmitt et al., and uploaded
# primary_cohort_TAD_boundaries.tgz to $WRKDIR/../ras_modifiers/data/Schmitt_2016/
# Also need to have pre-processed data (see main code for RAS modifiers)
while read gene; do
  # Upstream boundary
  zcat $WRKDIR/../ras_modifiers/data/Schmitt_2016/hg38/*.gz \
  | grep -E '^chr7' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools closest -D ref -io -id -a $WRKDIR/refs/genes.bed.gz -b - \
  | fgrep -w $gene | cut -f6 | sort -nk6,6 | head -n1
  # Downstream boundary
  zcat $WRKDIR/../ras_modifiers/data/Schmitt_2016/hg38/*.gz \
  | grep -E '^chr7' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools closest -D ref -io -iu -a $WRKDIR/refs/genes.bed.gz -b - \
  | fgrep -w $gene | cut -f6 | sort -nk7,7 | tail -n1
done < $WRKDIR/refs/genes.list | paste - -


### Find haplotype block decay around genes
mkdir $WRKDIR/data/1000G/
while read contig; do
  for suffix in gz gz.tbi; do
    cat << EOF > $WRKDIR/UGER/download_1000G_$suffix.sh
#!/usr/bin/env bash
wget -P $WRKDIR/data/1000G/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.$contig.filtered.SNV_INDEL_SV_phased_panel.vcf.$suffix
EOF
    chmod a+x $WRKDIR/UGER/download_1000G_$suffix.sh
    qsub -l h_rt=01:00:00 -V -N "download_1000G_$contig.$suffix" -cwd \
      $WRKDIR/UGER/download_1000G_$suffix.sh
    done
done < <( zcat $WRKDIR/refs/genes.bed.gz | cut -f1 )
# Define regions ±5Mb to query
bedtools slop \
  -b 5000000 \
  -i $WRKDIR/refs/genes.bed.gz \
  -g $WRKDIR/../ras_modifiers/refs/hg38.genome \
| bgzip -c > $WRKDIR/refs/genes.5Mb_buffer.bed.gz
# Process pedigrees and split into populations
wget -P $WRKDIR/data/1000G/ \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
for pop in AFR AMR EAS EUR SAS; do
  awk -v pop=$pop '{ if ($NF==pop && $3==0 && $4==0) print $2 }' \
    $WRKDIR/data/1000G/20130606_g1k_3202_samples_ped_population.txt \
  > $WRKDIR/data/1000G/$pop.samples.list
done
sed '1d' $WRKDIR/data/1000G/20130606_g1k_3202_samples_ped_population.txt \
| awk -v OFS="\t" '{ print $1, $2, $3, $4, $5, "-9" }' \
> $WRKDIR/data/1000G/1000G.ped
# Write job array to subset, sort, and index VCFs for each population then compute LD
for pop in AFR AMR EAS EUR SAS; do
  while read contig; do
    echo -e "$pop\t$contig"
  done < <( zcat $WRKDIR/refs/genes.bed.gz | cut -f1 )
done > $WRKDIR/UGER/1000G_LD.inputs.tsv
cat << EOF > $WRKDIR/UGER/1000G_LD.sh
#!/usr/bin/env bash

pop=\$1
contig=\$2

echo -e "Running contig \$contig for population \$pop"

export BASEDIR=/broad/VanAllenLab/xchip/cga_home/rcollins
export WRKDIR=/broad/VanAllenLab/xchip/cga_home/rcollins/egfr_modifiers
export CODEDIR=$WRKDIR/code/egfr_modifiers

source /broad/software/scripts/useuse
source /home/unix/rcollins/.bashrc

# Subset, sort, and index VCFs
time tabix -h \
  --regions $WRKDIR/refs/genes.5Mb_buffer.bed.gz \
  $WRKDIR/data/1000G/1kGP_high_coverage_Illumina.\$contig.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
| bcftools view \
  -O z -o $WRKDIR/data/1000G/1000G_phase3_highCov.\$pop.\$contig.RAS_loci.vcf.gz \
  --force-samples -S $WRKDIR/data/1000G/\$pop.samples.list \
  --min-ac 2 --max-alleles 2
tabix -f -p vcf $WRKDIR/data/1000G/1000G_phase3_highCov.\$pop.\$contig.RAS_loci.vcf.gz

# Compute LD
time plink \
  --r2 gz yes-really --ld-window-kb 5000 --ld-window 99999 --ld-window-r2 0.2 \
  --memory 4000 \
  --vcf $WRKDIR/data/1000G/1000G_phase3_highCov.\$pop.\$contig.RAS_loci.vcf.gz \
  --out $WRKDIR/data/1000G/1000G_phase3_highCov.\$pop.\$contig.RAS_loci
mv $WRKDIR/data/1000G/1000G_phase3_highCov.\$pop.\$contig.RAS_loci.ld.gz \
$WRKDIR/data/1000G/1000G_phase3_highCov.\$pop.\$contig.RAS_loci.ld.gzip \
&& zcat $WRKDIR/data/1000G/1000G_phase3_highCov.\$pop.\$contig.RAS_loci.ld.gzip \
| awk -v OFS="\t" '{ print \$1, \$2, \$3, \$4, \$5, \$6, \$7 }' | bgzip -c -l 5 \
> $WRKDIR/data/1000G/1000G_phase3_highCov.\$pop.\$contig.RAS_loci.ld.gz
tabix -f --sequence 1 --begin 2 --end 5 --skip-lines 1 \
  $WRKDIR/data/1000G/1000G_phase3_highCov.\$pop.\$contig.RAS_loci.ld.gz
EOF
chmod a+x $WRKDIR/UGER/1000G_LD.sh
cat << EOF > $WRKDIR/UGER/submit_1000G_LD.batch.sh
#!/usr/bin/env bash
#
#$ -t 1-5
#$ -l h_vmem=4G
#$ -l h_rt=01:00:00
#$ -o /broad/VanAllenLab/xchip/cga_home/rcollins/egfr_modifiers/UGER/logs/
#$ -e /broad/VanAllenLab/xchip/cga_home/rcollins/egfr_modifiers/UGER/logs/

export WRKDIR=/broad/VanAllenLab/xchip/cga_home/rcollins/egfr_modifiers

pop=\$( sed -n "\${SGE_TASK_ID}p" $WRKDIR/UGER/1000G_LD.inputs.tsv | cut -f1 )
contig=\$( sed -n "\${SGE_TASK_ID}p" $WRKDIR/UGER/1000G_LD.inputs.tsv | cut -f2 )

$WRKDIR/UGER/1000G_LD.sh \$pop \$contig
EOF
chmod a+x $WRKDIR/UGER/submit_1000G_LD.batch.sh
qsub $WRKDIR/UGER/submit_1000G_LD.batch.sh
# Define subset of variants with population-specific MAF >= 5%
for pop in AFR AMR EAS EUR SAS; do
  while read contig; do
    bcftools query \
      -i "INFO/AF_${pop} >= 0.05" \
      -f "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%AF_${pop}\n" \
      $WRKDIR/data/1000G/1000G_phase3_highCov.$pop.$contig.RAS_loci.vcf.gz
  done < <( zcat $WRKDIR/refs/genes.bed.gz | cut -f1 ) \
  | gzip -c > $WRKDIR/data/1000G/$pop.common_variants.tsv.gz
done
# Extract LD for all SNP-pairs with at least one SNP inside one of the three genes of interest
while read contig start end gene; do
  for pop in AFR AMR EAS EUR SAS; do
    zcat $WRKDIR/data/1000G/$pop.common_variants.tsv.gz \
    | awk -v contig=$contig -v start=$start -v end=$end \
      '{ if ($1==contig && $2>=start && $2<=end) print $5 }' \
    | fgrep -wf - <( zcat $WRKDIR/data/1000G/1000G_phase3_highCov.$pop.$contig.RAS_loci.ld.gz ) \
    | awk -v OFS="\n" '{ if ($7>=0.8) print $2, $5 }'
  done | sort -n | sed -e 1b -e '$!d' | paste -s -
done < <( zcat $WRKDIR/refs/genes.bed.gz )


