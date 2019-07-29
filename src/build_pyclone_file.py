import pandas as pd
import numpy as np
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('sample_name',help = 'sample name')
parser.add_argument('union_symbols',type = argparse.FileType('r'),help = 'path to list of SNVs in tsv format annotated with gene symbols')
parser.add_argument('samtools_snv',type = argparse.FileType('r'),help = 'path to output file from convert_pileup.py in tsv format')
parser.add_argument('titan_file',type = argparse.FileType('r'),help='path to segs.txt output file from TitanCNA')
args = parser.parse_args()



name = args.sample_name
union = pd.read_csv(args.union_symbols,sep='\t',names = ['SNV','SYMBOL'])
samtools_snv = pd.read_csv(args.samtools_snv.name, sep ='\t')
titan_file = pd.read_csv(args.titan_file.name,sep ='\t')




pyclonefile = samtools_snv.merge(union, on=['SNV'])



pyclonefile['mutation_id'] = pyclonefile['SYMBOL']+ ':' + pyclonefile['SNV'].str.split(':').str[0] +':'+ pyclonefile['SNV'].str.split(':').str[1]
pyclonefile = pyclonefile[['mutation_id','DP','ALT_COUNT','VAF']]
pyclonefile.columns = ['mutation_id','DP','var_counts','variant_freq']

pyclonefile['ref_counts'] = pyclonefile['DP'] - pyclonefile['var_counts']
pyclonefile = pyclonefile[['mutation_id','ref_counts','var_counts','variant_freq']]

pyclonefile['chrom'] = pyclonefile['mutation_id'].str.split(':').str[1]
pyclonefile['pos'] = pyclonefile['mutation_id'].str.split(':').str[2]
pyclonefile['normal_cn'] = 2
pyclonefile['minor_cn'] = 1
pyclonefile['major_cn'] = 1


pyclonefile.loc[pyclonefile['chrom'] == 'X','normal_cn'] = 1
pyclonefile.loc[pyclonefile['chrom'] == 'X','minor_cn'] = 0
pyclonefile['pos'] = pd.to_numeric(pyclonefile['pos'])



for chrom in range(1,22):
    titan_sub = titan_file[titan_file['Chromosome'] == chrom ]
    pyclone_sub = pyclonefile[pyclonefile['chrom'] == str(chrom)]
    start = titan_sub['Start_Position.bp.'].tolist()
    stop = titan_sub['End_Position.bp.'].tolist()
    major = titan_sub['MajorCN'].tolist()
    minor = titan_sub['MinorCN'].tolist()
    criteria = [pyclone_sub['pos'].between(i,j) for i , j in zip(start,stop)]


    pyclonefile.loc[:,'major_cn'][pyclonefile['chrom'] == str(chrom)] = np.select(criteria, major, 1)
    pyclonefile.loc[:,'minor_cn'][pyclonefile['chrom'] == str(chrom)] = np.select(criteria, minor,1)



pyclonefile = pyclonefile[['mutation_id','ref_counts','var_counts','normal_cn', 'minor_cn','major_cn','variant_freq']]

pyclonefile.to_csv(name + '_pyclone.tsv',sep='\t')
