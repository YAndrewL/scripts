#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 12:30:04 2019

@author: Yufan A. Liu
"""
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
#import numpy as np
import matplotlib.colors as colors
import numpy as np



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


#%%



#样本1
wt_1 = pd.read_table("N2-L1-1.anti.txt",header = None)
wt_1 = wt_1.iloc[:-5,:]
wt_1.columns = ['gene','counts']

wt_1['bin'] = 0
wt_1['rpkm'] = 0
wt_1_sum = np.sum(wt_1['counts'])
for i in range(len(wt_1)):
    wt_1['bin'][i] = list(gtf[gtf['gene'].isin([wt_1.gene[i]])].bin)[0]
    wt_1['rpkm'][i] = round((wt_1['counts'][i] * 1e9) / (list(gtf[gtf['gene'].isin([wt_1.gene[i]])].length)[0] * wt_1_sum),3)



#样本2 
wt_2 = pd.read_table("N2-L1-2.anti.txt",header = None)
wt_2 = wt_2.iloc[:-5,:]
wt_2.columns = ['gene','counts']
wt_2['bin'] = 0
wt_2['rpkm'] = 0
wt_2_sum = np.sum(wt_2['counts'])
for i in range(len(wt_2)):
    wt_2['bin'][i] = list(gtf[gtf['gene'].isin([wt_2.gene[i]])].bin)[0]
    wt_2['rpkm'][i] = round((wt_2['counts'][i] * 1e9) / (list(gtf[gtf['gene'].isin([wt_2.gene[i]])].length)[0] * wt_2_sum),3)
 
#样本3
    
    
    
cas501_1 = pd.read_table("cas501-1.anti.txt",header = None)
cas501_1 = cas501_1.iloc[:-5,:]
cas501_1.columns = ['gene','counts']
cas501_1['bin'] = 0
cas501_1['rpkm'] = 0
cas501_1_sum = np.sum(cas501_1['counts'])
for i in range(len(cas501_1)):
    cas501_1['bin'][i] = list(gtf[gtf['gene'].isin([cas501_1.gene[i]])].bin)[0]
    cas501_1['rpkm'][i] = round((cas501_1['counts'][i] * 1e9) / (list(gtf[gtf['gene'].isin([cas501_1.gene[i]])].length)[0] * cas501_1_sum),3)
    
#样本4
cas501_2 = pd.read_table("cas501-2.anti.txt",header = None)
cas501_2 = cas501_2.iloc[:-5,:]
cas501_2.columns = ['gene','counts']
cas501_2['bin'] = 0
cas501_2['rpkm'] = 0
cas501_2_sum = np.sum(cas501_2['counts'])
for i in range(len(cas501_1)):
    cas501_2['bin'][i] = list(gtf[gtf['gene'].isin([cas501_2.gene[i]])].bin)[0]
    cas501_2['rpkm'][i] = round((cas501_2['counts'][i] * 1e9) / (list(gtf[gtf['gene'].isin([cas501_2.gene[i]])].length)[0] * cas501_2_sum),3)
    


#分类汇总
'''
#按照rpkm
wt1_c = wt_1.groupby(by = ['bin']).agg({'rpkm':'sum'})
wt2_c = wt_2.groupby(by = ['bin']).agg({'rpkm':'sum'})
cas1_c = cas501_1.groupby(by = ['bin']).agg({'rpkm':'sum'})
cas2_c = cas501_2.groupby(by = ['bin']).agg({'rpkm':'sum'})
'''
wt1_c = wt_1.groupby(by = ['bin']).agg({'counts':'sum'})
wt2_c = wt_2.groupby(by = ['bin']).agg({'counts':'sum'})
cas1_c = cas501_1.groupby(by = ['bin']).agg({'counts':'sum'})
cas2_c = cas501_2.groupby(by = ['bin']).agg({'counts':'sum'})


heat_map = pd.merge(left = wt1_c, right = cas2_c, left_index = True, right_index = True)
heat_map.columns = ['N2','cas501']
heat_map['N2'] = np.log10(heat_map['N2'])
heat_map['cas501'] = np.log10(heat_map['cas501'])


#fig,ax = plt.subplots(figsize = (50,50))
ax = sns.clustermap(heat_map.T,xticklabels = True,figsize = (20,20)).savefig("antisense.pdf")

               #cmap =colors.LinearSegmentedColormap.from_list('Mycmap', ['#2D3E66','#ECEEE0','#5F4081'])).savefig("test.pdf")
               #cmap =colors.LinearSegmentedColormap.from_list('Mycmap', ['#FF3264','#6464FF']))
#ax.set_ytickslabel('')
#fig.tight_layout()
#2D3E66
#fig.savefig("test.pdf")
    
    
cheat = heat_map.drop('III-8')
ax = sns.clustermap(cheat.T,xticklabels = True,z_score = 0,figsize = (20,20))
#.savefig("cheat.pdf") 


    

    
    
    
    
    
    
    
    
    
    
    
    