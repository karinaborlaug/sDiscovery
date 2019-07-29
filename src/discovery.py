#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Main functinality for discovery of SNVs."""

import argparse
import parameters
import frames
import filter
import os
import pandas as pd


def main():
    """ Program containing the main logic to generate a set of tsv files.

    Args:
        param1: config_file

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', help='config_file')
    args = parser.parse_args()
    parameters.get_result_dir(args.config_file)
    strict_filter = parameters.get_strict_filter_parameters(args.config_file)
    loose_filter = parameters.get_loose_filter_parameters(args.config_file)
    vcf_files = parameters.get_vcf_files(args.config_file)
    maf_files = parameters.get_maf_files(args.config_file)
    types = parameters.get_type(args.config_file)
    path = 'results'

    # create sample objects
    samples = frames.make_sample_instances(vcf_files, maf_files)

    # remove tmp files
    for fname in os.listdir():
        if "tmp" in fname:
            os.remove(fname)

    # make results directory
    if not os.path.exists(path):
        os.mkdir(path)

    # create filters defined by user input for each sample
    for sample in samples:

        sample.set_strict_filter(strict_filter)
        sample.set_loose_filter(loose_filter)
        strict_list = sample.create_strict_filtered_snv()
        loose_list = sample.create_loose_filtered_snv()

    union_strict = filter.get_union_strict_filtered(samples)
    union_loose = filter.get_union_loose_filtered(samples)

    false_positives_strict = filter.get_false_positives(samples, union_strict, strict_filter['minimum_vaf_tumor'])
    false_positives_loose = filter.get_false_positives(samples, union_loose, loose_filter['minimum_vaf_tumor'])

    # set false positives for each sample
    for sample in samples:

        sample.set_false_positives_strict_filtered(false_positives_strict)
        sample.set_false_positives_loose_filtered(false_positives_loose)
        sample.set_strict_filter_snv_no_false_positives()
        sample.get_strict_filter_snv_no_false_positives().to_csv(os.path.join(path, sample.get_name() + "_strict.tsv"), sep='\t', index=False)
        sample.get_loose_filter_snv_no_false_positives().to_csv(os.path.join(path, sample.get_name() + "_loose.tsv"), sep='\t', index=False)

    union_strict_no_fp = union_strict[(~union_strict.SNV.isin(false_positives_strict.SNV))]

    union_loose_no_fp = union_loose[(~union_loose.SNV.isin(false_positives_loose.SNV))]

    common_loose_and_strict = union_loose_no_fp.merge(union_strict_no_fp, on=['SNV'])
    union_loose_no_fp_unique = union_loose_no_fp[(~union_loose_no_fp.SNV.isin(common_loose_and_strict.SNV))]
    union = pd.concat([union_strict_no_fp, union_loose_no_fp_unique], ignore_index=True)

    false_positives_strict.to_csv(os.path.join(path, 'false_positives_strict.tsv'), sep='\t', index=False, header=False)
    false_positives_loose.to_csv(os.path.join(path, 'false_positives_loose.tsv'), sep='\t', index=False)

    union_strict_no_fp[['SNV', 'SYMBOL']].to_csv(os.path.join('results', 'union_strict_with_symbol.tsv'), sep='\t', index=False, header=False)
    union_strict_no_fp[['SNV']].to_csv(os.path.join(path, 'union_strict.tsv'), sep='\t', index=False, header=False)

    union[['SNV', 'SYMBOL']].to_csv(os.path.join(path, 'union_loose_with_symbol.tsv'), sep='\t', index=False, header=False)
    union[['SNV']].to_csv(os.path.join(path, 'union_loose.tsv'), sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()
