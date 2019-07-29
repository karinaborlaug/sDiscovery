#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Filter functions for sample instances."""

import pandas as pd


def filter_snv(format_field, filter):
    """Test the format_field from SNV on a filter.

    Args:
        format_field:format_field from a vcf file
        filter: dictionary containing filtering thresholds

    Returns:
        The return value. True for success, False otherwise.

    """
    passed = False
    splitted = format_field.split(':')
    GT = splitted[0]
    AD = splitted[1]
    DP = int(splitted[3])

    AD_splitted = AD.split(",")
    alternate_read = float(AD_splitted[1])

    if GT == '0/1':
        vaf = alternate_read * 100 / DP
        if DP >= filter['minimum_tumor_coverage'] and alternate_read >= filter['minimum_alt_read_tumor'] and vaf >= filter['minimum_vaf_tumor']:
            passed = True
    else:
        if DP >= filter['minimum_normal_coverage'] and alternate_read <= filter['maximum_alt_read_normal']:
            passed = True

    return passed


def filter_ind(format_field, tum_cov, nor_cov, alt_read_tum, alt_read_nor, vaf_th, type):
    """Test the format_field from indel on a filter.

    Args:
        tum_cov: minimum tumor coverage
        nor_cov: minimum normal coverage
        alt_read_tum: minimum alternative read tumor
        alt_read_nor: maximum alternative read normal
        vaf_th: minimum vaf
        type: type of sample

    Returns:
        The return value. True for success, False otherwise.

    """
    passed = False
    splitted = format_field.split(':')
    DP = int(splitted[0])
    DP2 = splitted[1]
    TAR = splitted[2]
    TIR = splitted[3]
    TOR = splitted[4]
    DP50 = splitted[5]
    FDP50 = splitted[6]
    SUBDP50 = splitted[7]

    alt_read = int(TIR.split(',')[0])

    ref_read = DP - alt_read

    if type == 'tumor':

        vaf = alt_read * 100 / DP
        if DP >= 15 and alt_read >= 5 and vaf >= 5:
            passed = True
    if type == 'normal':

        if DP >= 10 and alt_read <= 1:
            passed = True

    return passed


def not_exonic(variant_class):
    """Check if SNV is exonic.

    Args:
        variant_class: the type of mutation

    Returns:
        True for not exonic, False otherwise.

    """
    variant_classes = ["5'Flank", 'Intron', 'RNA', "3'Flank", "3'UTR", "5'UTR", 'IGR']

    if variant_class in variant_classes:
        return False
    else:
        return True


def get_union_strict_filtered(samples):
    """Create union of strict filtered SNV.

    Args:
        samples: list of sample instances.

    Returns:
        union.

    """
    strict_filtered_snv = []
    for sample in samples:
        strict_filtered_snv.append(sample.get_strict_filter_snv())

    strict_filtered_snv_combined = pd.concat(strict_filtered_snv)

    union = strict_filtered_snv_combined[['SNV', 'SYMBOL']].drop_duplicates(subset='SNV', keep='first', inplace=False)
    union = union.reset_index(drop=True)
    return union


def get_union_loose_filtered(samples):
    """Create union of loose filtered SNV.

    Args:
        samples: list of sample instances.

    Returns:
        union.

    """
    loose_filtered_snv = []
    for sample in samples:
        loose_filtered_snv.append(sample.get_loose_filter_snv())

    loose_filtered_snv_combined = pd.concat(loose_filtered_snv)

    union = loose_filtered_snv_combined[['SNV', 'SYMBOL']].drop_duplicates(subset='SNV', keep='first', inplace=False)
    union = union.reset_index(drop=True)
    return union


def is_rejected_with_high_vaf(vaf, status, vaf_th):
    """Check if SNV is rejected with vaf over threshold.

    Args:
        vaf: vaf of the SNV.
        vaf_th: vaf threshold.
        status: reject/pass status of SNV.

    Returns:
        True if rejected, False otherwise

    """
    if vaf > vaf_th and status == 'REJECT':
        return True
    else:
        return False


def get_false_positives(samples, union, vaf_th):
    """Check if SNV is rejected with vaf over threshold.

    Args:
        samples: list of sample instances.
        union: unio of SNVs
        vaf_th: vaf threshold.


    Returns:
        True if rejected, False otherwise

    """
    for sample in samples:
        union[sample.get_name()] = False
        for i in range(len(union) - 1):
            if sample.snv_is_in_snv_list(union.iloc[i, 0]):
                vaf = sample.get_snv_vaf(union.iloc[i, 0])
                status = sample.get_snv_status(union.iloc[i, 0])
                false_positive = is_rejected_with_high_vaf(vaf, status, vaf_th)
                column_index = union.columns.get_loc(sample.get_name())
                union.iloc[i, column_index] = false_positive

    false_positives = pd.DataFrame()
    for col in union.columns[2:len(union.columns)]:
        false_positives = false_positives.append(union[union[col]])

    false_positives = false_positives.drop_duplicates(subset='SNV', keep='first', inplace=False)

    return false_positives


if __name__ == '__main__':
    d = {'minimum_tumor_coverage': 15, 'minimum_normal_coverage': 10, 'maksimum_alt_read_normal': 1, 'minimum_alt_read_tumor': 5, 'minimum_vaf_tumor': 5}
    print(filter_snv("0/1:109,20:80:70", d))
