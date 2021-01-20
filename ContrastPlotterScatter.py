# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd

bed_file = "gene.fil.bed"
f1 = "sta.txt"
f2 = "no_sta.txt"
gene_selected = "cdk-2"



locus_file = pd.read_table(bed_file,sep = " ", header = None)
gene_name = list(locus_file.iloc[:,10])
start_site = locus_file.iloc[:,1]
stop_site = locus_file.iloc[:,2]
# cdk-2的散点图

   
   #处理输入文件列表，使用1-base的数据   
header = ['Starve', 'No_starve']
f1_table = pd.read_table(f1,sep = "\t")
f1_table['ratio'] = round(f1_table['AD'] / f1_table['DP'] , 2 )
tb1 = f1_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)
f2_table = pd.read_table(f2,sep = "\t")
f2_table['ratio'] = round(f2_table['AD'] / f2_table['DP'] , 2 )
tb2 = f2_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis = 1)

   
pos = gene_name.index(gene_selected)
chrom = 'chr' + str(locus_file.iloc[pos,0])
start = start_site[pos]
stop = stop_site[pos]

tb1_gene = tb1[(tb1['#Chrom'] == chrom) & (tb1['Start(0base)'] >= start) & (tb1['Start(0base)'] <= stop)].drop(['#Chrom'], axis = 1)
tb2_gene = tb2[(tb2['#Chrom'] == chrom) & (tb2['Start(0base)'] >= start) & (tb2['Start(0base)'] <= stop)].drop(['#Chrom'], axis = 1)

#获取绘图列表并进行索引填充
tb1_gene = tb1_gene.set_index(['Start(0base)'])
tb2_gene = tb2_gene.set_index(['Start(0base)'])

for i in range(start,stop):
    if i not in tb1_gene.index:
        tb1_gene.loc[i] = -1
    if i not in tb2_gene.index:
        tb2_gene.loc[i] = -1
tb1_gene.sort_index(inplace = True) 
tb2_gene.sort_index(inplace = True)
   
   #合并
heat_table = pd.merge(tb1_gene,tb2_gene,left_index = True,right_index = True)
heat_table.columns = header
heat_table = heat_table.T

#%%
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
t1 = ax.scatter(heat_table.columns, heat_table.loc['Starve'],c = '#FC42B6',s = 20)
t2 = ax.scatter(heat_table.columns, heat_table.loc['No_starve'],c = '#5086FF',s = 20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel('Editing Level')
plt.ylim([0,1.1])
ax.legend((t1,t2),('Starve','No_starve'),loc = 'best')
plt.savefig("CDK2_scatter.pdf")








