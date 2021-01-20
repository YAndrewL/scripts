#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 13:57:49 2019

@author: jingyinghuang
"""

import pandas as pd
import os
import numpy as np
import seaborn as sns

os.chdir('/Users/jingyinghuang/Seafile/huangjingying/data2019/daf2/RNAediting/filter_res/')
##--------testing
#file='WT-1_RES_res.filter.txt'
#sample=file[:9]
#df=pd.read_csv(file,sep="\t",header=1)
#df.columns
#
##--get ratio for each site
#df['res_reads']=pd.DataFrame(df['AD:DP'].str.split(":",expand=True))[0].astype(int)
#df['all_reads']=pd.DataFrame(df['AD:DP'].str.split(":",expand=True))[1].astype(int)
#df['res_ratio']=df['res_reads']/df['all_reads']
#
##--get sum ratio of each 1Mbp
#df['length']=1000000
#df['ranges']=df['Start(0base)']//df['length']
#df['ranges'].drop_duplicates()
#df['chr_ranges']=df['#Chrom']+'-'+df['ranges'].astype(str)
#
#df_sum=pd.pivot_table(df,index=["chr_ranges"],values='res_ratio',aggfunc=['sum'])
#df_sum.columns=[sample]
#df_sum_merge=pd.merge(df_sum,df_sum,left_index=True,right_index=True,how='outer')#on='chr_ranges'

plot_path="../heatmap/"
name='GOU304'
name='GOU404'
name='GOU305'
files= os.listdir("./") 
files1=[i for i in files if name in i]
files2=[i for i in files if 'WT' in i]
files=files1+files2
print(files)
i=0
for i in range(len(files)):
    file=files[i]
    print(file)
    sample=file[:9]
    df=pd.read_csv(file,sep="\t",header=1)
    print(df.columns)
    
    if ('AD:DP' in df.columns):
        print('AD:DP')
        #--get ratio for each site
        df['res_reads']=pd.DataFrame(df['AD:DP'].str.split(":",expand=True))[0].astype(int)
        df['all_reads']=pd.DataFrame(df['AD:DP'].str.split(":",expand=True))[1].astype(int)
        df['res_ratio']=df['res_reads']/df['all_reads']
    elif ('AD' in df.columns and 'DP' in df.columns):
        print('AD','BD')
        df['res_ratio']=df['AD']/df['DP']
    #--get sum ratio of each 1Mbp
    df['length']=1000000
    df['ranges']=df['Start(0base)']//df['length']
    df['ranges'].drop_duplicates()
    df['chr_ranges']=df['#Chrom']+'-'+df['ranges'].astype(str)
    
    df_sum=pd.pivot_table(df,index=["chr_ranges"],values='res_ratio',aggfunc=['count'])
    df_sum.columns=[sample]

    if i==0:
        df_sum_merge=df_sum
    else:
        df_sum_merge=pd.merge(df_sum_merge,df_sum,left_index=True,right_index=True,how='outer')#on='chr_ranges'

df_sum_merge.fillna(0,inplace=True)
#df_sum_merge2=df_sum_merge/1000
df_sum_merge2=np.log10(df_sum_merge+1)
#--plot heatmap
df_sum_merge2.columns
#df_sum_merge2=df_sum_merge2[[ 'GOU3041-2', 'GOU3041-3','GOU3041-4', 'GOU3042-2','GOU3042-3','GOU3042-4', 
#        'WT-2_RES_', 'WT-1_RES_', 'WT-4_RES_']]
hier_allgenes = sns.clustermap(df_sum_merge2,col_cluster=True,row_cluster=True,
                               figsize=(8,20),
                                #cmap='seismic',
                                cmap="YlGnBu",
                                linewidths=.1,linecolor='grey',
#                                xticklabels =['WT','15min','30min','60min','90min'],
#                                vmax=1, vmin=0,
                                #cbar_kws={"orientation": "horizontal"},
                                yticklabels=df_sum_merge2.index,
                                )
hier_allgenes.savefig(plot_path+"heatmap_resCounts_"+str(name)+".pdf",format='pdf')
