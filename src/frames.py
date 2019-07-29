#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Converts vcf and maf files into sample instances."""

from sample import Sample
import pandas as pd
import io


def remove_comments_from_file(infile, outfile):
    """Remove header from VCF anf MAF files and writes to file.

    Args:
        infile: file to read.
        outfile: file to write.
    """
    with open(infile, 'r') as inf, open(outfile, 'w') as out:
        if infile.endswith('vcf'):
            for line in inf:
                if not line.startswith('##'):
                    out.write(line)
        else:
            for line in inf:
                if not line.startswith('#'):
                    out.write(line)


def file_to_pandas_dataframe(infile):
    """Convert file to pandas DataFrame.

    Args:
        infile: file to read.

    Returns:
            pandas DataFrame

    """
    with open(infile, 'r') as inf:
        if infile.endswith('vcf'):

            lines = [l for l in inf if not l.startswith('##')]
            dataframe = pd.read_csv(io.StringIO(''.join(lines)), sep='\t', dtype=str).rename(columns={'#CHROM': 'CHROM'})
        else:
            dataframe = pd.read_csv(inf, sep='\t', comment='#', dtype=str)

    return dataframe


def make_sample_instances(vcf_files, maf_files):
    """Make sample instances.

    Args:
        vcf_files: list of vcf files.
        maf_files: liest of maf files

    Returns:
            sample instances

    """
    samples = []
    for sample_id in vcf_files:
        remove_comments_from_file(vcf_files[sample_id], sample_id + "_tmp.vcf")
        remove_comments_from_file(maf_files[sample_id], sample_id + "_tmp.maf")
        vcf_df = file_to_pandas_dataframe(sample_id + "_tmp.vcf")
        maf_df = file_to_pandas_dataframe(sample_id + "_tmp.maf")
        vcf_df['SNV'] = vcf_df['CHROM'].astype(str) + ':' + vcf_df['POS'].astype(str) + ':' + vcf_df['REF'].astype(str) + ':' + vcf_df['ALT'].astype(str)
        snv = pd.concat([vcf_df, maf_df], axis=1)
        snv = snv.iloc[:, [0, 1, 2, 3, 4, 6, 9, 10, 11, 20, 51, 52, 53, 54, 55, 56, 72]]
        snv = snv[~(snv['CHROM'].astype(str).str.startswith('G'))]
        snv.name = sample_id
        samples.append(Sample(snv))
    return samples
