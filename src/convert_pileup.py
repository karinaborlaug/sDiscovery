import sys
import os
import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
import upsetplot
import seaborn as sns


sns.set()



parser = argparse.ArgumentParser()
parser.add_argument('union', type = argparse.FileType('r'),help = 'path to list of SNVs in tsv format annotated with gene symbols')
parser.add_argument('ref_pileup',nargs = 1,type = argparse.FileType('r'),help = 'file with samtools mpilup output in txt format.txt for the comparison sample ')
parser.add_argument('pileup_files',nargs = '+',type = argparse.FileType('r'),help='file(s) with samtools mpilup output in txt format,normal sample must be named normal.txt')
parser.add_argument('--min_alt',nargs = '?',type = int,default = 2,help = 'minimum alt_count,default = 2')
args = parser.parse_args()







def generate_result_frames(union,pileup,name,path):
    pileup = pd.read_csv(pileup, sep= '\t',names=['CHROM', 'POS','REF','DP','PILEUP','baseQ','mapQ','basePosOnReads'] ,header=None)
    pileup.name = name


    # Find out which ALT we are looking for
    samtools_result = pd.concat([pileup,union],ignore_index=False,axis=1, sort=False)

    samtools_result['ALT'] = samtools_result['SNV'].str.split(':').str[-1]

    #count the ALT_COUNT and calculate VAF


    samtools_result['ALT_COUNT'] = samtools_result.apply(lambda L : L.PILEUP.count(L.ALT) + L.PILEUP.count(L.ALT.lower()) if not pd.isnull(L.PILEUP) else 0, axis =1)
    samtools_result['VAF'] = samtools_result['ALT_COUNT'] / samtools_result['DP'] * 100
    samtools_result.reset_index(drop=True,inplace=True)

    #create samtools result
    samtools_result[['SNV','DP','ALT_COUNT','VAF']].to_csv(os.path.join(path, name +'_samtools_result.tsv'),sep = '\t',index = False)
    samtools_result = samtools_result[['SYMBOL','SNV','DP','ALT_COUNT','VAF']]



    samtools_result = samtools_result.set_index('SNV')


    return samtools_result

def generate_summary_stats(rapport,names,ref_name,min_alt,path):
    samples_dicts = {}
    for name in names:
        samples_dicts[name] = {}
        filtered_list = rapport[names]
        filtered_list = filtered_list.loc[rapport[name].str.split('/').str[0].astype(int) > min_alt]


        samples_dicts[name]['count'] = filtered_list[name].count()




        filtered_list.loc[:,'vaf-'+ name] = filtered_list[name].str.split('/').str[0].astype(int) / filtered_list[name].str.split('/').str[1].astype(int) * 100
        rapport.loc[:,'vaf-'+ name] = round(rapport[name].str.split('/').str[0].astype(int) / rapport[name].str.split('/').str[1].astype(int) * 100, 2)
        samples_dicts[name]['snv'] = filtered_list
        samples_dicts[name]['total_vaf']= float(filtered_list['vaf-'+ name].sum())
        samples_dicts[name]['mean'] = round(float(filtered_list['vaf-' + name].mean()),1)
        samples_dicts[name]['median']  = round(float(filtered_list['vaf-' + name].median()),1)
        samples_dicts[name]['std'] = round(float(filtered_list['vaf-' + name].std()),1)
        samples_dicts[name]['min']   = round(float(filtered_list['vaf-' + name].min()),1)
        samples_dicts[name]['max']  = round(float(filtered_list['vaf-' + name].max()),1)
        samples_dicts[name]['range'] = str(samples_dicts[name]['min'] ) + '-' + str(samples_dicts[name]['max'])



        names_to_use = samples_dicts[ref_name]['snv'].columns.difference(samples_dicts[name]['snv'].columns).tolist()
        merged = samples_dicts[ref_name]['snv'][names_to_use].merge(samples_dicts[name]['snv'],on=['SNV'])
        samples_dicts[name]['overlap'] = len(merged)
        if ((samples_dicts[name]['overlap'] == 0) or (samples_dicts[name]['count'] == 0)):
            samples_dicts[name]['percent'] = 0
        else: samples_dicts[name]['percent'] = round((samples_dicts[name]['overlap'] / samples_dicts[name]['count'] * 100),1)




    d1 = {'sample': names\
    ,'Total snv': [samples_dicts[names[i]]['count'] for i in range(len(names))]\
    ,'snv overlap with cfdna': [samples_dicts[names[i]]['overlap'] for i in range(len(names))]\
    ,'% overlap with cfdna': [samples_dicts[names[i]]['percent'] for i in range(len(names))]\
    ,'mean VAF%':[samples_dicts[names[i]]['mean'] for i in range(len(names))]\
    ,'median VAF%':[samples_dicts[names[i]]['median'] for i in range(len(names))]\
    ,'range VAF%':[samples_dicts[names[i]]['range'] for i in range(len(names))]\
    ,'std VAF%':[samples_dicts[names[i]]['std'] for i in range(len(names))]
    }
    df2 = pd.DataFrame(data=d1)



    ranks = df2.rank(axis=0, method='average', numeric_only=None, na_option='keep', ascending=False, pct=False)
    df2['rank snv overlap'] = ranks.iloc[:,2]
    df2['rank % overlap '] = ranks.iloc[:,3]

    rapport = rapport.set_index('GENE:POS')

    rapport.to_csv(os.path.join(path,'rapport.tsv'),sep ='\t',index = True)
    df2.to_csv(os.path.join(path,'stats.tsv'),sep ='\t',index = False)

    return df2



def generate_upsetplot(rapport,names,min_alt,path):
    # remove cases with no supporint reads in all samples from alt > min_alt
    rapport = rapport[rapport[names].apply(lambda alt : False if (alt.str.split('/').str[0].astype(int) < min_alt).all() else True,axis = 1)]


    upsetframe = rapport[names].reset_index()


    for name in names:


        mask  = (upsetframe[name].str.split('/').str[0].astype(int) >= min_alt)
        upsetframe.loc[mask,name] = True
        upsetframe.loc[upsetframe[name] != True,name] = False


    samples = [c for c in upsetframe.columns if c != 'SNV']
    samples_count_series = upsetframe.fillna(False).groupby(samples).count()['SNV']

    upsetplot.plot(samples_count_series,sort_by='cardinality')
    current_figure = plt.gcf()
    plt.title("Overlaps strict filtered SNVs",fontsize = 15)
    plt.ylabel("SNV count")

    current_figure.savefig(os.path.join(path,'upsetplot_' + str(min_alt) + '.png'))




def main():
    union = args.union
    ref = args.ref_pileup[0]
    pileups = args.pileup_files
    min_alt = args.min_alt

    #union = pd.read_csv(union,sep='\t',names=['SNV'],header = None)
    union = pd.read_csv(union,sep='\t',names=['SNV','SYMBOL'],header = None)

    #union = union['SNV']

    name_ref = os.path.abspath(ref.name)
    name_ref = os.path.basename(name_ref)
    name_ref = ref_name = name_ref.split('.')[0]
    path='rapports_' + str(min_alt)

    # make rapports directory
    if not os.path.exists(path):
        os.mkdir(path)


    ref_frame = generate_result_frames(union,ref,name_ref,path)
    ref_frame[name_ref] = ref_frame['ALT_COUNT'].astype(str) +'/'+ ref_frame['DP'].astype(str)


    frames = [ref_frame[name_ref]]
    names = [name_ref]
    for p in pileups:
        name = os.path.abspath(p.name)
        name = os.path.basename(name)
        name = name.split('.')[0]
        result_frame = generate_result_frames(union,p,name,path)
        result_frame['GENE:POS'] = result_frame['SYMBOL']+ '|' + result_frame.index.str.split(':').str[0] +':'+ result_frame.index.str.split(':').str[1]

        result_frame[name] = result_frame['ALT_COUNT'].astype(str) +'/'+ result_frame['DP'].astype(str)


        result_frame[name].head()
        frames.append(pd.concat([result_frame[name],result_frame['GENE:POS']],axis = 1))
        names.append(name)



    rapport = pd.concat(frames, axis=1)
    rapport = rapport.loc[:,~rapport.columns.duplicated()]



    #filter out snv with no supporting reads in any of the samples
    rapport = rapport[rapport[names].apply(lambda alt : False if (alt.str.split('/').str[0].astype(int) < 2).all() else True,axis = 1)]
    # filter out SNV with alt_read > 2 in normal if normal is provided
    if "normal" in names:
        rapport = rapport[rapport["normal"].apply(lambda alt : False if int(alt.split('/')[0]) > 1 else True)]







    generate_upsetplot(rapport,names,min_alt,path)
    summary_stats = generate_summary_stats(rapport,names,ref_name,min_alt,path)



main()
