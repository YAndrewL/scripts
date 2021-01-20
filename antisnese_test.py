#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:02:55 2020

@author: Yufan A. Liu
"""
import pandas as pd
import sys
import seaborn as sns
#import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

#调试用
bed_file = 'gene.fil.bed'
f1 = 'EX5048.txt'
f2 = 'GOU3806.txt'
gene_selected = 'dyf-5'
id1 = 'EX5048'
id2 = 'GOU3086'

#处理bed文件
locus_file = pd.read_table(bed_file,sep = " ", header = None)
gene_name = list(locus_file.iloc[:,10])
#start_site = locus_file.iloc[:,1]
#stop_site = locus_file.iloc[:,2]

#处理输入文件列表，使用1-base的数据   
header = [id1, id2]
f1_table = pd.read_table(f1,sep = "\t")

f1_table['ratio'] = round(f1_table['AD'] / f1_table['DP'] , 2 )
tb1 = f1_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'Unnamed: 8','AD', 'DP'], axis = 1)
f2_table = pd.read_table(f2,sep = "\t")
f2_table['ratio'] = round(f2_table['AD'] / f2_table['DP'] , 2 )
tb2 = f2_table.drop(['End(1base)', 'Supporting_reads', 'Strand', 'Type', 'AD','Unnamed: 8', 'DP'], axis = 1)

#提取gene所在的染色体位置

pos = gene_name.index(gene_selected)
chrom = 'chr' + str(locus_file.iloc[pos,0])
start = 9358845
stop = 9364210
#start = start_site[pos]
#stop = stop_site[pos]

tb1_gene = tb1[(tb1['#Chrom'] == chrom) & (tb1['Start(0base)'] >= start) & (tb1['Start(0base)'] <= stop)].drop(['#Chrom'], axis = 1)
tb2_gene = tb2[(tb2['#Chrom'] == chrom) & (tb2['Start(0base)'] >= start) & (tb2['Start(0base)'] <= stop)].drop(['#Chrom'], axis = 1)
 
#获取绘图列表并进行索引填充
tb1_gene = tb1_gene.set_index(['Start(0base)'])
tb2_gene = tb2_gene.set_index(['Start(0base)'])

for i in range(start,stop):
    if i not in tb1_gene.index:
        tb1_gene.loc[i] = 0
    if i not in tb2_gene.index:
        tb2_gene.loc[i] = 0
tb1_gene.sort_index(inplace = True) 
tb2_gene.sort_index(inplace = True)

#合并
heat_table = pd.merge(tb1_gene,tb2_gene,left_index = True,right_index = True)
heat_table.columns = header
heat_table = heat_table.T

#绘图
fig, ax = plt.subplots(figsize=(20, 2.5))
font = {'family' : 'Times New Roman', 'Weight' : 'normal', 'Size': 15,}
sns.heatmap(heat_table,cmap =colors.LinearSegmentedColormap.from_list('Mycmap', ['#A6A4A4', '#FC0006']), vmax = 0.6, vmin = 0)

ax.set_title(gene_selected + " editing ratio graph." + chrom + " positon:" + str(start) + ":" +str(stop), font)                                                                                 
ax.set_xlabel("Position", font)
ax.set_ylabel("Sample ID", font)
#ax.set_yticks(range(-0.5, 1.5))

  #  tick_locator = ticker.MaxNLocator(6)
  #  ax.xaxis.set_major_locator(tick_locator) #可以显示的最多的横轴坐标数目
#ax.set_yticklabels(header, va = "center",rotation = 45)
   # ax.set_xticklabels(rotation = 45)
#plt.xticks(range(start, stop, int((stop-start)/5)))
#ax.set_xticks(np.arange(start, stop, int((stop-start)/5)))
#ax.set_xticks([])
ax.axhline(1, linewidth=2, c='black')

fig.tight_layout()   


fig.savefig(gene_selected + ".pdf")