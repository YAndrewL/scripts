#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 15:05:00 2019

@author: Yufan A. Liu
"""

#Editing gene annotation script
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import scipy

#def main():
#读入gene.flt.bed文件提取基因名称

locus_file = pd.read_table("gene.fil.bed",sep = " ", header = None)
gene_name = list(locus_file.iloc[:,10])
start_site = locus_file.iloc[:,1]
stop_site = locus_file.iloc[:,2]

#读入带有生物学重复的样本数据,度量值是每个基因内的编辑数目
header = ['N2_1','N2_2','cas501_1','cas501_2' ]

f1_table = pd.read_table("N2_1.res",sep = "\t")
f1_table['ratio'] = round(f1_table['AD'] / f1_table['DP'] , 2 )
#tb1 = f1_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)

f2_table = pd.read_table("N2_2.res",sep = "\t")
f2_table['ratio'] = round(f2_table['AD'] / f2_table['DP'] , 2 )
#tb2 = f2_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)

f3_table = pd.read_table("cas501_1.res",sep = "\t")
f3_table['ratio'] = round(f3_table['AD'] / f3_table['DP'] , 2 )
#tb2 = f2_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)

f4_table = pd.read_table("cas501_2.res",sep = "\t")
f4_table['ratio'] = round(f4_table['AD'] / f4_table['DP'] , 2 )
#tb2 = f2_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)


#确定每个基因和基因范围

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
data = {'N2_1':s1, 'N2_2':s2, 'cas501_1':s3, 'cas501_2':s4} #使用字典通过Series生成Dataframe
sum_table = pd.DataFrame(data)
flt_table = sum_table[~(sum_table == 0).all(axis =1)] #选择在样本中编辑数目不全为0的点

'''
#评估样本间大概服从什么分布,方差远大于均值，近似使用负二项分布NB
cor_data = pd.DataFrame(columns = ['length','number','var','mean'],index = flt_table.index)
for gene in flt_table.index:
    pos = gene_name.index(gene)
    start = start_site[pos]
    stop = stop_site[pos]
    gene_length = int(stop - start)
    N2_num = flt_table['N2_1'][gene]
    cor_data['length'][gene] = gene_length
    cor_data['number'][gene] = N2_num
    cor_data['mean'][gene] = (flt_table['N2_1'][gene] + flt_table['N2_2'][gene]) / 2
    mean = (flt_table['N2_1'][gene] + flt_table['N2_2'][gene]) / 2
    cor_data['var'][gene] = (np.square((flt_table['N2_1'][gene] - mean)) + np.square((flt_table['N2_2'][gene] - mean)))/2
#sns.scatterplot(cor_data[cor_data['length']<200000]['length'],cor_data[cor_data['number']< 1200]['number'])
sns.scatterplot(cor_data['var'],cor_data['mean'])
'''

#flt_table['DiffVarNorm'] = (flt_table.iloc[:,2:4].mean(1) - flt_table.iloc[:,0:2].mean(1)) / flt_table.iloc[:, 0:4].var(1)
flt_table['DiffChang'] = flt_table.iloc[:,2:4].mean(1) - flt_table.iloc[:,0:2].mean(1)
#w, pval = scipy.stats.wilcoxon(flt_table.iloc[:,2:4].mean(1),flt_table.iloc[:,0:2].mean(1) )
#flt_name = list(flt_name)

#假设检验，对于已经编辑的样本，假设样本的数据来自于一个正态分布
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
flt_table.loc[flt_table['DiffChang'] >= 50,'sig'] = 'up'
flt_table.loc[flt_table['DiffChang'] <= -50,'sig'] = 'down'

p_flt_table = flt_table[flt_table['pvalue'] >=1]
p_flt_table['ChangeAdj'] = np.log10(p_flt_table['DiffChang'])

'''
fig,ax = plt.subplots()
font = {'Family' : 'Times New Roman'}
sns.scatterplot(x = 'DiffChang', y = 'pvalue',hue = 'sig', hue_order=('Down', 'Normal', 'Up'), palette=("blue", "gray", "red"), data = flt_table)
ax.set_xlabel("NormDiffChange",font)
ax.set_ylabel("-log(Pvalue)",font)
'''
#绘制差异表达的火山图

heat_table['sig'] = 'normal'
heat_table.loc[heat_table['DiffChang'] >= 50,'sig'] = 'up'
heat_table.loc[heat_table['DiffChang'] <= -50,'sig'] = 'down'
#heat_table = p_flt_table.drop(['N2_1','N2_2','cas501_1','cas501_2','ChangeAdj','sig','DiffVarNorm'],axis = 1)
fig ,ax = plt.subplots()
font = {'Family' : 'Times New Roman'}
sns.scatterplot('DiffChang', 'pvalue', data = heat_table,hue = 'sig', hue_order = ('down','normal','up'), palette=('blue', 'gray', 'red'),s=20)
ax.set_xlabel('Correlative Difference', font)
ax.set_ylabel('-log(pvalue)', font)
ax.set_title('Volcano plot',font)
ax.text(400,1.9,'dyf-5')
fig.tight_layout()
fig.savefig("VolcanoPlotofDyf-5.pdf")






#main()