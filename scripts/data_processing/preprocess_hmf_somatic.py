#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

"""
Combine HMF somatic SNVs/indels and CNAs and subset to gene(s) of interest
"""


import argparse
import pandas as pd
import pybedtools as pbt
import pysam
from re import sub
from sys import stdout


chrom_dict = {}
for contig in range(23):
    chrom_dict[str(contig)] = contig
    chrom_dict['chr{}'.format(contig)] = contig
chrom_dict.update({'X' : 22, 'chrX' : 22, 'Y' : 23, 'chrY' : 23})


def build_cnv_records(cna_in, gtf_in, vcf_out, key_genes=None):
    """
    Creates one VCF record for each gene CNA listed in cna_in
    """

    # Load all CNAs
    cna_df = pd.read_csv(cna_in, sep='\t', header=None)
    cna_df.columns = 'sample gene cna'.split()

    # Get list of genes to consider
    genes = set(cna_df.gene.values.tolist())
    if key_genes is not None:
        genes = genes.intersection(set(key_genes))
    if len(genes) > 0:
        genes = sorted(list(genes))
    else:
        msg = 'Failed to identify any gene symbols shared between --genes-gtf ' + \
              'and --priority-genes. Are you sure you are providing the correct ' + \
              'files as input?'
        exit(msg)

    # Filter GTF to genes of interest
    gbt = pbt.BedTool(gtf_in).filter(lambda f: f[2] == 'gene').\
                              filter(lambda f: f.attrs.get('gene_name') in genes).\
                              sort()

    # Emit one record per CNA per gene
    cna_records = []
    for gdat in gbt:

        gene = gdat.attrs.get('gene_name')
        g_cna = cna_df.loc[cna_df.gene == gene, :]
        for svtype in g_cna.cna.value_counts().index:

            # Make basic record
            chrom = sub('^chr', '', gdat.chrom)
            pos = gdat.start
            end = gdat.stop
            alleles = ('N', '<' + svtype + '>', )
            vid = '{}_{}_{}_{}'.format(chrom, pos, *alleles)
            record = vcf_out.new_record(contig=chrom, start=pos, stop=end,
                                        alleles=alleles, id=vid)

            # Populate samples
            AC, AN, AF = 0, 0, 0
            cna_samps = set(g_cna.loc[g_cna.cna == svtype, 'sample'].values.tolist())
            for sid in record.samples.keys():
                AN += 1
                if sid in cna_samps:
                    GT = 1
                    AC += 1
                else:
                    GT = 0
                record.samples[sid]['GT'] = GT

            # Update AC/AN/AF in INFO
            record.info['AN'] = AN
            record.info['AC'] = AC
            record.info['AF'] = AC / AN

            # Add record to list of records to add to output VCF
            cna_records.append(record)

    return cna_records


def clean_mut_record(old_record, vcf_out):
    """
    Cleans a VCF record for a single somatic mutation
    """

    # Get reused values from old record
    vid = '{}_{}_{}_{}'.format(old_record.chrom, old_record.pos, 
                               old_record.ref, old_record.alleles[1])
    star_alleles = [idx for idx, a in enumerate(old_record.alleles) if a == '*']
    new_alleles = (a for idx, a in enumerate(old_record.alleles) if a != '*')

    # Make new record attached to vcf_out header
    # Dev note: for reasons I can't fully explain, it appears that trying to
    # write a record from vcf_in (with a subtly different header) to vcf_out
    # fails due to a suspicious htslib error. This error goes away when 
    # vcf_out is opened in 'wb' mode, and the header entry cited by the error 
    # changes depending on input and output files. I suspect this is more 
    # htslib/pysam weirdness but have been unable to track it down.
    # Long story short: it seems easier to make a new record native to the 
    # header of vcf_out and transfer the info over from old_record
    record = vcf_out.new_record(contig=old_record.chrom, start=old_record.pos,
                                alleles=new_alleles, id=vid, qual=old_record.qual)
    if 'END' in record.info.keys():
        record.info.pop('END')
    record.filter.add('PASS')

    # Translate all GTs to allele dosages (presence|absence)
    AC, AN, AF = 0, 0, 0
    for sid, sdat in old_record.samples.items():
        GT = []
        for a in sdat['GT']:
            if a is None:
                continue
            elif a in star_alleles:
                GT.append(0)
            else:
                GT.append(a)
        if len(GT) == 0:
            new_GT = (None, )
        else:
            AN += 1
            if any([a > 0 for a in GT]):
                new_GT = (1, )
                AC += 1
            else:
                new_GT = (0, )
        record.samples[sid]['GT'] = new_GT

    # Update AC/AN/AF
    record.info['AC'] = AC
    record.info['AN'] = AN
    if AN > 0:
        record.info['AF'] = round(AC / AN, 8)
    else:
        record.info['AF'] = 0
    
    return record


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--muts-vcf', metavar='VCF', required=True,
                        help='Somatic mutation VCF. Must include all samples.')
    parser.add_argument('--cna-tsv', help='Three-column .tsv of sample ID, gene, ' +
                        'and CNA (DEL|AMP). If included, will be added to output VCF.')
    parser.add_argument('--genes-gtf', help='.gtf of gene annotations', required=True)
    parser.add_argument('--priority-genes', help='If provided, will subset all data ' + 
                        'to only variants in these genes [default: keep all genes]')
    parser.add_argument('--ref-fasta', help='Reference .fasta', required=True)
    parser.add_argument('--header', help='Header to use for output .vcf', required=True)
    parser.add_argument('-o', '--outfile', default='stdout', help='path to somatic ' +
                        'variation .vcf [default: stdout]')
    args = parser.parse_args()

    # Load genes of interest, if optioned
    if args.priority_genes is not None:
        with open(args.priority_genes) as fin:
            key_genes = list(set([l.rstrip() for l in fin.readlines()]))
    else:
        key_genes = None

    # Load VCF header
    header = pysam.VariantFile(args.header).header

    # Open connections to in/out VCFs
    vcf_in = pysam.VariantFile(args.muts_vcf)
    for sid in vcf_in.header.samples:
        header.add_sample(sid)
    vcf_out = pysam.VariantFile(args.outfile, 'w', header=header)

    # Build VCF records for all CNAs, if provided
    if args.cna_tsv is not None:
        cna_records = build_cnv_records(args.cna_tsv, args.genes_gtf, 
                                        vcf_out, key_genes)
    else:
        cna_records = []

    # Iterate over mutation records, spiking in CNA records where appropriate
    for record in vcf_in.fetch():
        mut_chrom = record.chrom
        mut_pos = record.pos

        # Check if any CNA records should be spiked into output VCF
        if len(cna_records) > 0:
            next_cna = cna_records[0]
            cna_chrom = next_cna.chrom
            cna_pos = next_cna.pos
            while chrom_dict[cna_chrom] <= chrom_dict[mut_chrom] \
            and cna_pos <= mut_pos:
                vcf_out.write(next_cna)

                # Move to next CNA record (if any exist)
                if len(cna_records) > 1:
                    cna_records = cna_records[1:]
                    next_cna = cna_records[0]
                    cna_chrom = next_cna.chrom
                    cna_pos = next_cna.pos
                else:
                    break

        # After handling CNA records, clean somatic mutation record 
        # and write to --outfile
        record = clean_mut_record(record, vcf_out)
        vcf_out.write(record)

    # After processing all records, close connection to output VCF to clear buffer
    vcf_out.close()


if __name__ == '__main__':
    main()

