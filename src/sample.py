#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""class sample to contain SNVs and filtering functions."""
import pandas as pd
import filter


class Sample:
    """Class to conatain SNV info and filtering functinality.

    Args:
        snv_list: pandas DataFrame with SNVs
        indels_list: pandas DataFrame with indels
        false_positives_strict_filter: pandas DataFrame false positives strict filtered
        false_positives_loose_filter: pandas DataFrame false positives loose filtee
        strict_filtered_snv: pandas DataFrame with strict filtered SNVs
        strict_filtered_snv_no_false_positives: pandas DataFrame strict filtered SNVs with no FP
        passed_exonic: pandas DataFrame with passed exonic SNV

    """

    def __init__(self,snv_list = None,indels_list = None):
        self.snv_list = snv_list
        self.snv_list['VAF'] = pd.to_numeric(self.snv_list['t_alt_count']) * 100 / pd.to_numeric(self.snv_list['t_depth'])
        self.indels_list = indels_list
        self.false_positives_strict_filter = 0
        self.false_positives_loose_filter = 0
        self.strict_filtered_snv = 0
        self.strict_filtered_snv_no_false_positives = 0
        self.passed_exonic = 0

    def set_strict_filter(self, strict_filter):
        self.strict_filter = strict_filter

    def set_loose_filter(self, loose_filter):
        self.loose_filter = loose_filter

    def set_false_positives_strict_filtered(self, false_positives):
        self.false_positives_strict_filter = false_positives
        self.set_strict_filter_snv_no_false_positives()

    def set_false_positives_loose_filtered(self, false_positives):
        self.false_positives_loose_filter = false_positives
        self.set_loose_filter_snv_no_false_positives()

    def set_strict_filter_snv_no_false_positives(self):
        self.strict_filtered_snv_no_false_positives = self.strict_filtered_snv[(~self.strict_filtered_snv.SNV.isin(self.false_positives_strict_filter.SNV))]

    def set_loose_filter_snv_no_false_positives(self):
        loose_filtered_snv_no_false_positives = self.loose_filtered_snv[(~self.loose_filtered_snv.SNV.isin(self.false_positives_loose_filter.SNV))]
        self.loose_filtered_snv_no_false_positives = pd.concat([loose_filtered_snv_no_false_positives, self.strict_filtered_snv_no_false_positives], ignore_index=True)

    def get_strict_filter_snv_no_false_positives(self):
        return self.strict_filtered_snv_no_false_positives

    def get_loose_filter_snv_no_false_positives(self):
        return self.loose_filtered_snv_no_false_positives

    def get_name(self):
        return self.snv_list.name

    def snv_is_in_snv_list(self, snv):
        return (self.snv_list['SNV'] == snv).any()

    def get_snv_list(self):
        return self.snv_list

    def get_snv_vaf(self, snv):
        x = self.snv_list[self.snv_list['SNV'] == snv].index.values.astype(int)[0]
        y = self.snv_list.columns.get_loc('VAF')
        return self.snv_list.iloc[x, y]

    def get_snv_status(self, snv):
        x = self.snv_list[self.snv_list['SNV'] == snv].index.values.astype(int)[0]
        y = self.snv_list.columns.get_loc('FILTER')
        return self.snv_list.iloc[x, y]

    def get_indels(self):
        return self.indels_list

    def get_snv_count(self):
        return self.snv_list.POS.count()

    def get_indel_count(self):
        return self.indels_list.POS.count()

    def get_strict_filter_snv(self):
        return self.strict_filtered_snv

    def get_loose_filter_snv(self):
        return self.loose_filtered_snv

    def create_strict_filtered_snv(self):
        passed = self.snv_list[self.snv_list['FILTER'] == 'PASS']
        passed_exonic = passed[passed['Variant_Classification'].apply(filter.not_exonic)]
        passed_exonic = passed_exonic.dropna(subset=['Variant_Classification'])

        filtered = passed_exonic[passed_exonic.iloc[:, 6].apply(filter.filter_snv, args=(self.strict_filter,))]
        filtered = filtered[filtered.iloc[:, 7].apply(filter.filter_snv, args=(self.strict_filter,))]

        self.strict_filtered_snv = filtered.loc[:, ['SNV', 't_depth', 't_alt_count', 'VAF', 'SYMBOL']]
        self.passed_exonic = passed_exonic

        return self.strict_filtered_snv

    def create_loose_filtered_snv(self):
        removed_strict_filtered = self.passed_exonic[(~self.passed_exonic.SNV.isin(self.strict_filtered_snv.SNV))]

        filtered = removed_strict_filtered[removed_strict_filtered .iloc[:, 6].apply(filter.filter_snv, args=(self.loose_filter,))]
        filtered = filtered[filtered.iloc[:, 7].apply(filter.filter_snv, args=(self.loose_filter,))]

        self.loose_filtered_snv = filtered.loc[:, ['SNV', 't_depth', 't_alt_count', 'VAF', 'SYMBOL']]

        return self.loose_filtered_snv

    def set_false_positives(self, false_positives):
        self.false_positives = false_positives
