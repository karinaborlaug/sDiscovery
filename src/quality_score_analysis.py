#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Program to analyse phred scores,create plots and list of thresholds."""

import sys
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker
from collections import Counter
import argparse
sns.set()


parser = argparse.ArgumentParser()
parser.add_argument('union', type=argparse.FileType('r'), help='path to list of SNVs in tsv format')
parser.add_argument('pileup_files', nargs='+', type=argparse.FileType('r'), help='file(s) with samtools mpilup output in txt format')
args = parser.parse_args()


def remove_chars(char):
    """removes unwanted chars from pileup string.

    Args:
        char: string from pileup

    Returns:
        filtered string

    """
    return re.sub('(\$|\^{1}.{1}|\*|\\+[0-9]+[ACGTNacgtn]+|-[0-9]+[ACGTNacgtn]+)', '', char)


# ref : https://stackoverflow.com/questions/11122291/how-to-find-char-in-string-and-get-all-the-indexes
def find_index_of_all_char_in_string(s, char):
    """find index of all chars in a string s.

    Args:
        char: ascii character
        s: string from fileup

    Returns:
        list of indexes

    """
    return [i for i, ltr in enumerate(s) if ltr.upper() == char]


def convert_phred_to_number(char):

    ASCII = { "!" : 33, "\"" : 34, "#" : 35, "$" : 36, "%" : 37 ,
            "&": 38, "\'": 39 ,"(" : 40 , ")": 41, "*" : 42, "+" :43,
            "," :44, "-":45 , "." :46, "/" : 47 ,"0" : 48,"1" : 49,
            "2": 50,"3" : 51,"4" :52,"5" :53,"6" : 54,"7" : 55,"8" : 56,
            "9" : 57, ":" : 58,";" : 59,"<" :60,"=" :61,">" :62,"?" :63,"@" : 64,
            "A" : 65,"B" : 66, "C" : 67, "D" : 68, "E" : 69, "F" : 70,
            "G" :71,"H" : 72,"I" : 73,"J":74,"K":75,"L" : 76 ,"M" : 77,
            "N" : 78,"O" : 79,"P": 80,"Q" : 81,"R" : 82,"S" : 83,"T" : 84,
            "U" :85 ,"V" : 86,"W" : 87,"X" : 88,"Y" :89,"Z" : 90, "[" : 91,
            "\\" : 92,"]" :93,"^" : 94,'_' : 95 ,"`" :96,"a" : 97,"b" : 98,
            "c" : 99,"d" : 100,"e" : 101, "f" : 102, "g" : 103 ,"h" : 104,
            "i" : 105, "j" : 106, "k" : 107,"l" : 108, "m" : 109, "n" :110,
            "o" : 111,"p" : 112, "q" : 113 ,"r" : 114 , "s" : 115,"t" : 116,
            "u" : 117 , "v" : 118 , "w" : 119 ,"x" : 120, "y" : 121, "z" :122,
            }

    return ASCII[char] - 33


# https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
def generate_frequency_dicts(result_frame):
    """generated frequency dictionarys of phred quality.

    Args:
        result_frame: pandas result frame from pileup

    Returns:
        baseQ frequency dictionar,mapQ frequency dictionary

    """
    baseQ = [result_frame.iloc[x, 13] for x in range(len(result_frame))]
    flat_list_baseQ = [item for baseQ in baseQ for item in baseQ]

    mapQ = [result_frame.iloc[x, 14] for x in range(len(result_frame))]
    flat_list_mapQ = [item for mapQ in mapQ for item in mapQ]

    baseQ_frequency_dict = {x: flat_list_baseQ.count(x) for x in flat_list_baseQ}
    mapQ_frequency_dict = {x: flat_list_mapQ.count(x) for x in flat_list_mapQ}

    return baseQ_frequency_dict, mapQ_frequency_dict


def generate_result_frames(pileup, union):
    """generate pandas result fram from pileup and snv union of sets.

    Args:
        pileup: samtools pileup
        union: list of SNVs

    Returns:
        result DataFrame

    """
    result_frame = pd.concat([pileup, union], ignore_index=False, axis=1, sort=False)
    result_frame['ALT'] = result_frame['SNV'].str.split(':').str[-1]
    result_frame['PILEUP'] = result_frame['PILEUP'].apply(remove_chars)
    result_frame['INDICES_ALT'] = result_frame.apply(lambda x: find_index_of_all_char_in_string(x['PILEUP'], x['ALT']), axis=1)
    result_frame.index.name = "mutation"
    result_frame['baseQ_ALT'] = result_frame.apply(lambda x: [x['baseQ'][i] for i in x['INDICES_ALT']], axis=1)
    result_frame['mapQ_ALT'] = result_frame.apply(lambda x: [x['mapQ'][i] for i in x['INDICES_ALT']], axis=1)
    result_frame['baseQ_ALT_score'] = result_frame.apply(lambda x: [convert_phred_to_number(x['baseQ_ALT'][i]) for i in range(len(x['baseQ_ALT']))], axis=1)
    result_frame['mapQ_ALT_score'] = result_frame.apply(lambda x: [convert_phred_to_number(x['mapQ_ALT'][i]) for i in range(len(x['mapQ_ALT']))], axis=1)

    return result_frame


def quality_analysis_plot(baseQ_frequency_dict, mapQ_frequency_dict, result_frame, name):
    """generate frequency and swarmplot.

    Args:
        baseQ_frequency_dict: frequency dictionary baseQ score
        mapQ_frequency_dict: frequency dictionary mapQ score
        result_frame: pandas DataFrame
        name: title for plot

    """

    # baseq frequency plot
    plt.figure()
    plt.bar(baseQ_frequency_dict.keys(), baseQ_frequency_dict.values(), color='g')
    plt.title(name + ': baseQ alternative reads', fontsize=15)
    plt.xlabel('baseQ-value', fontsize=15)
    plt.ylabel('SNV count', fontsize=15)
    plt.yscale('log',basey = 2)
    plt.savefig(os.path.join('q_scores', name + '_baseQ_freq.png'))

    # mapq frequency plot
    plt.figure()
    plt.bar(mapQ_frequency_dict.keys(), mapQ_frequency_dict.values(), color='g')
    plt.title(name + ': mapQ alternative reads', fontsize=15)
    plt.xlabel('mapQ-value', fontsize=15)
    plt.ylabel('SNV count', fontsize=15)
    plt.yscale('log',basey = 2)
    plt.savefig(os.path.join('q_scores', name + '_mapQ_freq.png'))

    # baseq per mutatation
    plt.figure()
    ax = sns.swarmplot(data=result_frame['baseQ_ALT_score'])
    plt.title(name + ': baseQ per mutation,alternative reads', fontsize=15)
    plt.xlabel('SNV', fontsize=15)
    plt.ylabel('baseQ-value', fontsize=15)

    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 5))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    plt.savefig(os.path.join('q_scores', name + '_baseQ_per_mut.png'))

    # mapq per mutatation
    plt.figure()
    ax = sns.swarmplot(data=result_frame['mapQ_ALT_score'])
    plt.title(name + ': mapQ per mutation,alternative reads', fontsize=15)
    plt.xlabel('SNV', fontsize=15)
    plt.ylabel('mapQ-value', fontsize=15)

    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, 5))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    plt.savefig(os.path.join('q_scores', name + '_mapQ_per_mut.png'))

    plt.close('all')


def generate_all_sample_dict(list_of_dicts):
    """generate a dictionary for with phred scores for all the samples.

    Args:
        list_of_dicts: list of dictionary per sample

    Returns:
        all samples dictionary

    """
    list_of_counters = []
    for dict in list_of_dicts:
        list_of_counters.append(Counter(dict))

    current_count = list_of_counters[0]
    for i in range(1, len(list_of_counters) - 1):
        current_count = current_count + list_of_counters[i]

    return current_count


def main():
    union = args.union
    pileups = args.pileup_files

    if not os.path.exists('q_scores'):
        os.mkdir('q_scores')

    union = pd.read_csv(union, sep='\t', names=['SNV'], header=None)

    list_of_baseQ_freq_dict = []
    list_of_mapQ_freq_dict = []
    total_frame = pd.DataFrame()

    names = []
    for p in pileups:
        name = os.path.abspath(p.name)
        name = os.path.basename(name)
        name = name.split('.')[0]
        names.append(name)
        pileup = pd.read_csv(p, sep='\t', names=['CHROM', 'POS', 'REF', 'DP', 'PILEUP', 'baseQ', 'mapQ', 'basePosOnReads'], header=None)

        result_frame = generate_result_frames(pileup, union)

        total_frame[name + '_baseQ'] = result_frame['baseQ_ALT_score']
        total_frame[name + '_mapQ'] = result_frame['mapQ_ALT_score']

        baseQ_frequency_dict, mapQ_frequency_dict = generate_frequency_dicts(result_frame)

        quality_analysis_plot(baseQ_frequency_dict, mapQ_frequency_dict, result_frame, name)
        list_of_baseQ_freq_dict.append(baseQ_frequency_dict)
        list_of_mapQ_freq_dict.append(mapQ_frequency_dict)

    total_frame['baseQ_ALT_score'] = np.empty((len(total_frame), 0)).tolist()
    total_frame['mapQ_ALT_score'] = np.empty((len(total_frame), 0)).tolist()

    for n in names:

        total_frame['baseQ_ALT_score'] = total_frame['baseQ_ALT_score'] + total_frame[n + '_baseQ']
        total_frame['mapQ_ALT_score'] = total_frame['mapQ_ALT_score'] + total_frame[n + '_mapQ']

    baseQ_frequency_dict_all = generate_all_sample_dict(list_of_baseQ_freq_dict)
    mapQ_frequency_dict_all = generate_all_sample_dict(list_of_mapQ_freq_dict)




    total_frame['baseQ_th'] = total_frame['baseQ_ALT_score'].apply(lambda x: 60 if not x else (5 if max(x) >= 30 else 20))
    total_frame['mapQ_th'] = total_frame['mapQ_ALT_score'].apply(lambda x: 60 if not x else (5 if max(x) >= 30 else 20))

    pd.Series(total_frame['baseQ_th']).to_csv(os.path.join('q_scores', "baseQ_values.tsv"), index=False, sep='\t', header=False)
    pd.Series(total_frame['mapQ_th']).to_csv(os.path.join('q_scores', "mapQ_values.tsv"), index=False, sep='\t', header=False)
    quality_analysis_plot(baseQ_frequency_dict_all, mapQ_frequency_dict_all, total_frame, "across_all_samples")


main()
