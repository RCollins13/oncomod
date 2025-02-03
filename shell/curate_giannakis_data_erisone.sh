#!/usr/bin/env bash

###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Curate data from Gurjao et al. 2021 for FGFR4:KRAS replication analysis

# Note: intended to be executed on the MGB ERISOne cluster


### Set local parameters
export BASEDIR=/data/gusev/PROFILE
export WRKDIR=/data/gusev/USERS/rlc47/giannakis
if [ ! -e $WRKDIR ]; then
  mkdir $WRKDIR
fi
export CODEDIR=$WRKDIR/../code/oncomod
export VEP_CACHE=$WRKDIR/../refs/vep_cache
export VEP_PLUGINS=$WRKDIR/../code/vep_plugins
export KGDIR=/data/gusev/1000G/
cd $WRKDIR


### Set up directory trees as necessary
for SUBDIR in data data/tumor_mafs data/normal_vcfs data/common_SNPs \
              data/sample_info \
              LSF LSF/scripts LSF/logs; do
  if ! [ -e $WRKDIR/$SUBDIR ]; then
    mkdir $WRKDIR/$SUBDIR
  fi
done


### Parse all tumor MAFs and only keep KRAS mutations
# Download & compress MAFs
while read sid; do
  awk -v sid="$sid" -v FS="\t" '{ if ($1==sid) print $10 }' \
    $WRKDIR/data/sample_info/giannakis_april_2023.pair.tsv
done < $WRKDIR/data/sample_info/gurjao_2021.900_published_crc_patients.pair_ids \
| gsutil -m cp -I $WRKDIR/data/tumor_mafs/
find $WRKDIR/data/tumor_mafs/ -name "*.maf" | xargs -I {} gzip {}
# Extract KRAS mutations from all MAFs
zcat $WRKDIR/data/tumor_mafs/*.maf.gz \
| awk -v FS="\t" -v OFS="\t" \
  '{ if ($1=="KRAS") print $16, $1, $36, $288, $5, $6, $11, $13, $256 }' \
| sort -Vk1,1 -k2,2V -k3,3V -k4,4V -k5,5V -k6,6V -k7,7V -k8,8V -k9,9V \
| cat <( echo -e "pair_id\tgene\tHUGO_protein\taa_csq\tchrom\tpos\tref\talt\tHGVS" ) - \
| gzip -c \
> $WRKDIR/data/gurjao_2021.kras_mutations.tsv.gz


### Curate germline FGFR4 variant data
# Download, compress, and index all DeepVariant VCFs
while read sid; do
  awk -v sid="$sid" -v FS="\t" '{ if ($1==sid) print $10 }' \
    $WRKDIR/data/sample_info/giannakis_april_2023.sample.tsv
done < <( fgrep Normal $WRKDIR/data/sample_info/gurjao_2021.900_published_crc_patients.exome_sample_ids.list ) \
| gsutil -m cp -I $WRKDIR/data/normal_vcfs/
find $WRKDIR/data/normal_vcfs/ -name "*.vcf" | xargs -I {} bgzip {}
find $WRKDIR/data/normal_vcfs/ -name "*.vcf.gz" | xargs -I {} tabix {}
# Merge single-sample VCFs into cohort-wide VCF for FGFR4 locus
find $WRKDIR/data/normal_vcfs/ -name "*.vcf.gz" \
> $WRKDIR/data/normal_vcfs/gurjao_2021.sample_vcfs.list
bcftools merge \
  --file-list $WRKDIR/data/normal_vcfs/gurjao_2021.sample_vcfs.list \
  --regions 5:176013887-177025145 \
  --no-version \
  --missing-to-ref \
  --threads 4 \
| bcftools annotate \
  -x FORMAT/AD,FORMAT/DP,FORMAT/PL \
| bcftools norm \
  --check-ref x \
  -m - \
  --fasta-ref /data/gusev/USERS/rlc47/TCGA/refs/GRCh37.fa \
| bcftools +fill-tags - -- -t AN,AC,AF \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
| bcftools view \
  --trim-alt-alleles \
  -Oz -o $WRKDIR/data/gurjao_2021.FGFR4.vcf.gz \
  --include 'AC > 0'
tabix -p vcf -f $WRKDIR/data/gurjao_2021.FGFR4.vcf.gz
# Annotate with VEP
/data/gusev/USERS/rlc47/code/ensembl-vep/vep \
  --input_file $WRKDIR/data/gurjao_2021.FGFR4.vcf.gz \
  --format vcf \
  --output_file $WRKDIR/data/gurjao_2021.FGFR4.vep.vcf \
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
  --fasta /data/gusev/USERS/rlc47/TCGA/refs/GRCh37.fa \
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
  --custom $VEP_CACHE/GTEx_eQTLs.vcf.gz,GTEx,vcf,exact,0,GTEx_eGene,GTEx_eQTL_beta,GTEx_eQTL_tissue
bgzip -f $WRKDIR/data/gurjao_2021.FGFR4.vep.vcf
tabix -p vcf $WRKDIR/data/gurjao_2021.FGFR4.vep.vcf.gz
# Clean VEP output
$CODEDIR/scripts/data_processing/cleanup_vep.py \
  --gtf $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
  --transcript-info $WRKDIR/../refs/gencode.v19.annotation.transcript_info.tsv.gz \
  --mode RAS_loci \
  --priority-genes <( echo "FGFR4" ) \
  $WRKDIR/data/gurjao_2021.FGFR4.vep.vcf.gz \
  stdout \
| grep -ve "^##bcftools" | grep -ve "^##CADD_" | grep -ve "^##UTRAnnotator_" \
| grep -ve "^##LoF_" | grep -ve "^##SpliceAI_" | grep -ve "^##VEP-command-line" \
| bgzip -c \
> $WRKDIR/data/gurjao_2021.FGFR4.vep.clean.vcf.gz
tabix -p vcf $WRKDIR/data/gurjao_2021.FGFR4.vep.clean.vcf.gz
# Apply blanket GQ<15 mask and only retain sites with at least one GT above this threshold
bcftools +setGT $WRKDIR/data/gurjao_2021.FGFR4.vep.clean.vcf.gz -- \
  -t q -n . -i 'GQ < 10' \
| bcftools +fill-tags -- -t AN,AC,AF \
| bcftools view -i 'INFO/AC > 0' \
  -Oz -o $WRKDIR/data/gurjao_2021.FGFR4.vep.clean.gq10.vcf.gz
tabix -p vcf -f $WRKDIR/data/gurjao_2021.FGFR4.vep.clean.gq10.vcf.gz
# Build allele dosage matrix
$CODEDIR/scripts/data_processing/vcf2dosage.py \
  $WRKDIR/data/gurjao_2021.FGFR4.vep.clean.gq10.vcf.gz - \
| gzip -c \
> $WRKDIR/data/gurjao_2021.FGFR4.vep.clean.gq10.dosage.tsv.gz
# Make table of codon position, protein change, and in-cohort frequency for all variants
echo -e "#vid\tcsq\thgsvp\taa_number\taa_change\tAF\tgnomAD_popmax" \
> $WRKDIR/data/gurjao_2021.FGFR4.gq10.variant_csqs.tsv
bcftools +split-vep \
  -f '%ID\t%SYMBOL\t%Consequence\t%HGVSp\t%Protein_position\t%Amino_acids\t%AF\t%INFO/gnomAD_AF_POPMAX\n' \
  $WRKDIR/data/gurjao_2021.FGFR4.vep.clean.gq10.vcf.gz \
| awk -v FS="\t" -v OFS="\t" '{ if ($2=="FGFR4") print }' | cut -f1,3- \
>> $WRKDIR/data/gurjao_2021.FGFR4.gq10.variant_csqs.tsv


### Compute PCs for normal samples
# Subset common 1000G SNPs to (very rough) coding regions likely captured on exomes
zcat $WRKDIR/../refs/gencode.v19.annotation.gtf.gz \
| awk -v OFS="\t" '{ if ($3=="exon") print $1, $4, $5 }' \
| sed 's/^chr//g' | sort -Vk1,1 -k2,2n -k3,3n | grep -e '^[1-9]' \
| bedtools merge -i - \
| bgzip -c \
> $WRKDIR/../refs/gencode.v19.exonic_regions.b37.bed.gz
bcftools view \
  --regions-file $WRKDIR/../refs/gencode.v19.exonic_regions.b37.bed.gz \
  -Oz -o $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.coding.pruned.vcf.gz \
  $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.pruned.vcf.gz
tabix -p vcf -f $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.coding.pruned.vcf.gz
for contig in $( seq 1 22 ); do
  bcftools view \
    --regions $contig \
    -O z -o $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.coding.pruned.$contig.vcf.gz \
    $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.coding.pruned.vcf.gz
  tabix -p vcf -f $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.coding.pruned.$contig.vcf.gz
done
# Find common coding SNPs that overlap with common 1000G SNPs
for contig in $( seq 1 22 ); do
  cat << EOF > $WRKDIR/LSF/scripts/extract_common_SNPs.$contig.sh
#!/usr/bin/env bash
. /PHShome/rlc47/.bashrc
cd $WRKDIR
bcftools merge \
  --file-list $WRKDIR/data/normal_vcfs/gurjao_2021.sample_vcfs.list \
  --regions-file $KGDIR/phase3/common_SNPs/LD_pruned/1000G.ph3.common_snps.coding.pruned.$contig.vcf.gz \
  --regions-overlap 2 \
  --no-version \
  --missing-to-ref \
  --threads 4 \
| bcftools +setGT -- -t q -n . -i 'GQ < 5' \
| bcftools annotate \
  -x FORMAT/AD,FORMAT/DP,FORMAT/PL,FORMAT/PL,FORMAT/GQ,FORMAT/VAF \
| bcftools norm \
  --check-ref x \
  -m - \
  --fasta-ref /data/gusev/USERS/rlc47/TCGA/refs/GRCh37.fa \
| bcftools +fill-tags - -- -t AN,AC,AF,F_MISSING \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
| bcftools view \
  --trim-alt-alleles \
  -Oz -o $WRKDIR/data/common_SNPs/gurjao_2021.common_SNPs.$contig.vcf.gz \
  --include 'INFO/AF > 0.01 & INFO/F_MISSING < 0.05'
tabix -p vcf -f $WRKDIR/data/common_SNPs/gurjao_2021.common_SNPs.$contig.vcf.gz
EOF
  chmod a+x $WRKDIR/LSF/scripts/extract_common_SNPs.$contig.sh
  rm $WRKDIR/LSF/logs/extract_common_SNPs.$contig.*
  bsub \
    -q normal -sla miket_sc \
    -n 4 -R 'rusage[mem=12000]' -J extract_common_SNPs_$contig \
    -o $WRKDIR/LSF/logs/extract_common_SNPs.$contig.log \
    -e $WRKDIR/LSF/logs/extract_common_SNPs.$contig.err \
    $WRKDIR/LSF/scripts/extract_common_SNPs.$contig.sh
done
# Merge common SNP data
for contig in $( seq 1 22 ); do
  echo $WRKDIR/data/common_SNPs/gurjao_2021.common_SNPs.$contig.vcf.gz
done > $WRKDIR/data/common_SNPs/gurjao_2021.common_SNPs.shards.list
bcftools concat -a \
  --file-list $WRKDIR/data/common_SNPs/gurjao_2021.common_SNPs.shards.list \
| bcftools norm \
  --atomize \
  --check-ref x \
  --fasta-ref /data/gusev/USERS/rlc47/TCGA/refs/GRCh37.fa \
  --threads 4 \
| bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
  -O z -o $WRKDIR/data/common_SNPs/gurjao_2021.common_SNPs.vcf.gz
tabix -p vcf -f $WRKDIR/data/common_SNPs/gurjao_2021.common_SNPs.vcf.gz
# Compute PCs from common SNPs
module load plink/2.0a2.3
plink2 \
  --threads 8 \
  --memory 32000 \
  --vcf $WRKDIR/data/common_SNPs/gurjao_2021.common_SNPs.vcf.gz \
  --geno 0.1 \
  --pca \
  --out $WRKDIR/data/gurjao_2021.common_SNPs



