#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 13:39:51 2019

@author: Yufan A. Liu
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import scipy

locus_file = pd.read_table("gene.fil.bed",sep = " ", header = None)
gene_name = list(locus_file.iloc[:,10])
start_site = locus_file.iloc[:,1]
stop_site = locus_file.iloc[:,2]

#读入带有生物学重复的样本数据,度量值是每个基因内的编辑数目
header = ['WT_1','WT_2','GOU3041_1','GOU3041_2' ]

f1_table = pd.read_table("WT-1_RES_res.filter.txt",sep = "\t")
f1_table2 = f1_table['AD:DP'].str.split(":",expand = True)
f1_table2.columns = ['AD','DP']
f1_table = f1_table.join(f1_table2).drop('AD:DP',axis =1)
#f1_table['ratio'] = round(f1_table['AD'].astype(int) / f1_table['DP'].astype(int) , 2 )
#tb1 = f1_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)

f2_table = pd.read_table("WT-2_RES_res.filter.txt",sep = "\t")
#f2_table['ratio'] = round(f2_table['AD'] / f2_table['DP'] , 2 )
f2_table2 = f2_table['AD:DP'].str.split(":",expand = True)
f2_table2.columns = ['AD','DP']
f2_table = f2_table.join(f2_table2).drop('AD:DP',axis =1)
#f2_table['ratio'] = round(f2_table['AD'] / f2_table['DP'] , 2 )
#tb2 = f2_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)

f3_table = pd.read_table("GOU3041-2_RES_res.filter.txt",sep = "\t")
f3_table['ratio'] = round(f3_table['AD'] / f3_table['DP'] , 2 )
#tb2 = f2_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)

f4_table = pd.read_table("GOU3041-3_RES_res.filter.txt",sep = "\t")
f4_table['ratio'] = round(f4_table['AD'] / f4_table['DP'] , 2 )
#tb2 = f2_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)

s1 = pd.Series()
s2 = pd.Series()
s3 = pd.Series()
s4 = pd.Series()
for gene in gene_name:
    pos = gene_name.index(gene)
    chrom = 'chr' + str(locus_file.iloc[pos,0])
    start  = start_site[pos]
    stop = stop_site[pos]
#    sum_num = len(f1_table[(f1_table['#Chrom'] == chrom) & (f1_table['Start(0base)'] >= start) & (f1_table['Start(0base)'] <= stop)])
#    s1[gene] = sum_num
#确定每个基因在每个样本中的编辑数目
    for i in range(4):
        tab = locals()["f" + str(i+1)+"_table"]
        sum_num = len(tab[(tab['#Chrom'] == chrom) & (tab['Start(0base)'] >= start) & (tab['Start(0base)'] <= stop)])
        locals()["s" + str(i+1)][gene] = sum_num                
data = {'WT_1':s1, 'WT_2':s2, 'GOU3041_1':s3, 'GOU3042_2':s4} #使用字典通过Series生成Dataframe
sum_table = pd.DataFrame(data)
flt_table = sum_table[~(sum_table == 0).all(axis =1)] #选择在样本中编辑数目不全为0的点
flt_table['DiffChang'] = flt_table.iloc[:,2:4].mean(1) - flt_table.iloc[:,0:2].mean(1)





pval_list = []
for i in range(len(flt_table.index)):
    w,pval = scipy.stats.ttest_ind(flt_table.iloc[i,2:4],flt_table.iloc[i,0:2])
    #,correction = True,zero_method = 'zsplit')
    if pval == 0 or pval == 'nan':
        pval =1
    p10 = -math.log10(pval)
    pval_list.append(p10)
flt_table['pvalue'] = pval_list
flt_table['sig'] = 'normal'
flt_table.loc[flt_table['DiffChang'] >= 20,'sig'] = 'up'
flt_table.loc[flt_table['DiffChang'] <= -20,'sig'] = 'down'

p_flt_table = flt_table[flt_table['pvalue'] >=1]
#p_flt_table['ChangeAdj'] = np.log10(p_flt_table['DiffChang'])

fig ,ax = plt.subplots()
font = {'Family' : 'Times New Roman'}
sns.scatterplot('DiffChang', 'pvalue', data = p_flt_table,hue = 'sig', hue_order = ('down','normal','up'), palette=('blue', 'gray', 'red'),s=20)
ax.set_xlabel('Correlative Difference', font)
ax.set_ylabel('-log(pvalue)', font)
ax.set_title('Volcano plot',font)
fig.tight_layout()
fig.savefig("GOU3041-WT.pdf")

p_flt_table.to_excel("GOU3041-WT.xls")


f2_table['AD'] = f2_table['AD'].astype('int')
f2_table['DP'] = f2_table['DP'].astype('int')
f2_table.to_csv("WT.txt",sep = "\t" ,index = 0)