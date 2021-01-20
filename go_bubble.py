#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:55:20 2020

@author: Yufan A. Liu
"""

import numpy as np  
import matplotlib.pyplot as plt  
import os
import pandas as pd
import seaborn as sns

#%% data
#os.chdir(r"G:\huangjingying\data\starve\get_TEresult\go")
go=pd.read_excel("go_cas501_up.xlsx",sheet_name="Sheet1" )
go.columns
#['Sublist', 'Category', 'Term', 'RT', 'Genes', 'Count', '%', 'P-Value',
#       'Benjamini']
#%% plot
cm = plt.cm.get_cmap('GnBu')
plt.figure(figsize=(4,8))
ax = plt.subplot(111)
test=plt.scatter(go['sublist'],go["Description"]
    ,s=go['%']*30,c=(-go['Log10(P)']),marker='o',alpha=1,cmap=cm,label=go['%']
    ,edgecolors='black')
#position=fig.add_axes([0.15, 0.05, 0.7, 0.03])#位置[左,下,右,上]
#plt.colorbar(cax=position,orientation='horizontal')
#plt.legend(loc=2, bbox_to_anchor=(1.2,0.85), ncol=1, frameon=True, fontsize=12,
#handlelength=2,  borderpad = 1.8,
#handletextpad=1, title='enrichment', scatterpoints = 1)

plt.colorbar(orientation='horizontal')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
plt.tick_params(
    axis='y',  # changes apply to the x-axis,axis : {‘x’, ‘y’, ‘both’} Axis on which to operate; default is ‘both’.
    which='both',  # both major and minor ticks are affected
    #bottom=False,  # ticks along the bottom edge are off
    #top=False,  # ticks along the top edge are off
    #labelbottom=False, # labels along the bottom edge are off
    left=False
    #labelleft=False
    )
#plt.xticks(['f30_vs_starve', 'f1h_vs_starve'],['f30_vs_starve', 'f1h_vs_starve'],fontsize=15)

plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
#plt.title(title,fontsize=15)
#plt.legend(loc=2,bbox_to_anchor=(1.2,0.85),borderaxespad = 0.,frameon=False,ncol=2,title='enrichment(%)',fontsize=15)
plt.savefig("go.png",dpi = 600, bbox_inches = 'tight')


#%%

test2=plt.scatter(go['sublist'],go["Description"]
    ,s=go['%']*30,c='white',marker='o',alpha=1,label=go['%']
    ,edgecolors='black')
plt.legend(loc=2, bbox_to_anchor=(1.2,0.85),borderaxespad = 0.,frameon=False,ncol=2,title='enrichment(%)',fontsize=15)
