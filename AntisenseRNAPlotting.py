#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 11:02:05 2019

@author: Yufan A. Liu
"""
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
#import numpy as np
import matplotlib.colors as colors
import numpy as np

#%% cell1
gtf = pd.read_table("Celegans.gtf",sep = '\t',header = None)
gtf.columns = ['chrom', 'site', 'gene','length']


gtf['gene'] = gtf['gene'].drop_duplicates()
gtf = gtf.fillna(0)

gtf = gtf[~gtf['gene'].isin([0])]
gtf.reset_index(drop = True, inplace = True)

bin_list = []
for i in range(len(gtf)):
    rang = gtf['site'][i] // 1000000
    bin_name = gtf['chrom'][i] + '-' + str(rang +1)
    bin_list.append(bin_name)
gtf['bin'] = bin_list
#%% cell2

count_table = pd.read_table("Norm.txt",sep = '\t',header = 0,index_col = 0)

#%%
tablebin = []
for i in range(len(count_table)):
    tablebin.append(list(gtf[gtf['gene'].isin([count_table.index[i]])].bin)[0])
count_table['bin'] = tablebin
#%% cell3

heat = count_table.groupby('bin').agg({'N2_1':['sum'], 'N2_2':['sum'],'N2_3':['sum'], 'cas501_1':['sum'],'cas501_2':['sum'],'cas501_3':['sum'],'GOU2348_1':['sum'],'GOU2348_2':['sum'],'GOU2348_3':['sum'],'GOU2387_1':['sum'],'GOU2387_2':['sum'],'GOU2387_3':['sum']})
heat.columns =['N2_1', 'N2_2','N2_3',  'cas501_1', 'cas501_2', 'cas501_3','GOU2348_1', 'GOU2348_2','GOU2348_3','GOU2387_1', 'GOU2387_2','GOU2387_3']

#heat = heat.drop(['N2_3','cas501_3'],axis = 1)

heat['N2'] = np.mean([heat['N2_1'],heat['N2_2'],heat['N2_3']],axis = 0)
heat['cas501'] = np.mean([heat['cas501_1'],heat['cas501_2'],heat['cas501_3']],axis = 0)
heat['GOU2348'] = np.mean([heat['GOU2348_1'],heat['GOU2348_2'],heat['GOU2348_3']],axis = 0)
heat['GOU2387'] = np.mean([heat['GOU2387_1'],heat['GOU2387_2'],heat['GOU2387_3']],axis = 0)

heat_map = np.log2(heat.iloc[:,-4:])
heat_map = heat_map.drop(['MtDNA-1'])
#heat_map = heat_map
#heat_map.to_csv("test.csv")
#%% cell4 drawing
#fig,ax = plt.subplots(figsize=(40,20))

#plt.figure(figsize = (40,20))

g = sns.clustermap(heat_map,xticklabels = True,
               #standard_scale=1,
               col_cluster = False,
               #ax = ax,
               #cbar_kws = {ax:'ax'},
               #cmap = "mako",
               figsize = (10,20),
               linewidths=.1,
               #cmap =colors.LinearSegmentedColormap.from_list('Mycmap', ['#7B4F9D','#F3F8ED','#08407D','#A1D3B0']),
               #map =colors.LinearSegmentedColormap.from_list('Mycmap', ['#08407D','#A1D3B0','#F3F8ED','#7B4F9D']),
               #方案2
               cmap = colors.LinearSegmentedColormap.from_list('Mycmap', ['#F3F8ED','#A1D3B0','#08407D','#7B4F9D']),
               #配色方案3
               yticklabels = heat_map.index,
               cbar_kws={'ticks':[10,12,14,16]})
               #vmax = 15, vmin = 0)
              #cmap =colors.LinearSegmentedColormap.from_list('Mycmap', ['#2D3E66','#ECEEE0','#5F4081']))


#plt.tight_layout()
font = {'family' : 'Times New Roman', 'Weight' : 'normal', 'Size': 15,}
g.ax_heatmap.set_ylabel("Chorm bins",font)      
g.ax_heatmap.set_xlabel("Samples",font)
#g.ax_heatmap.legend(loc = 'center left',labels = 'Legend' )
#g.ax_row_dendrogram.legend(loc="upper right", bbox_to_anchor=(0.47, 0.8),ncol=1)
#g.tight_layout()
g.savefig("antisense.pdf")


















