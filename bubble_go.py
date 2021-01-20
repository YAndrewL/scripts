import numpy as np  
import matplotlib.pyplot as plt  
import os
import pandas as pd
 
import seaborn as sns


os.chdir(r"G:\huangjingying\data\starve\get_TEresult\go")
title="biological process"
go=pd.read_excel("go_mrnanotchange_rpfup.xlsx",sheet_name=title )
go.columns
#['Sublist', 'Category', 'Term', 'RT', 'Genes', 'Count', '%', 'P-Value',
#       'Benjamini']

cm = plt.cm.get_cmap('GnBu_r')
plt.figure(figsize=(4,8))
ax = plt.subplot(111)
test=plt.scatter(go['Sublist'],go["Term"]
    ,s=go['%']*10,c=np.log10(go['P-Value']),marker='o',alpha=1,cmap=cm,label=go['%']
    ,edgecolors='black')
#position=fig.add_axes([0.15, 0.05, 0.7, 0.03])#位置[左,下,右,上]
#plt.colorbar(cax=position,orientation='horizontal')
#plt.legend(loc=2, bbox_to_anchor=(1.2,0.85), ncol=1, frameon=True, fontsize=12,
#handlelength=2,  borderpad = 1.8,
#handletextpad=1, title='enrichment', scatterpoints = 1)

#plt.legend()

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
plt.title(title,fontsize=15)
path=r"G:\huangjingying\data\starve\plots\\"
plt.savefig(path+title+'_go_mrnanotchange_rpfup+'+title+'.png',dpi=600, bbox_inches='tight')

#####--plot legand
go['%']
plt.figure(figsize=(4,8))
ax = plt.subplot(111)
np.min(go['%'])
n=7
test2=plt.scatter(go.loc[n,'Sublist'],go.loc[n,"Term"]
    ,s=go.loc[n,'%']*10,c='white',marker='o',alpha=1,label=go.loc[n,'%']
    ,edgecolors='black')

np.max(go['%'])
n=15
test2=plt.scatter(go.loc[n,'Sublist'],go.loc[n,"Term"]
    ,s=go.loc[n,'%']*10,c='white',marker='o',alpha=1,label=go.loc[n,'%']
    ,edgecolors='black')
plt.legend(loc=2, bbox_to_anchor=(1.2,0.85),borderaxespad = 0.,frameon=False,ncol=2,title='enrichment(%)',fontsize=15)
plt.savefig(path+title+'_go_mrnanotchange_rpfup+'+title+'_legand.png',dpi=600, bbox_inches='tight')

###################

title="KEGG pathway"
go=pd.read_excel("go_mrnanotchange_rpfup.xlsx",sheet_name=title )
go.columns
#['Sublist', 'Category', 'Term', 'RT', 'Genes', 'Count', '%', 'P-Value',
#       'Benjamini']

cm = plt.cm.get_cmap('GnBu_r')
plt.figure(figsize=(4,8))
ax = plt.subplot(111)
test=plt.scatter(go['Sublist'],go["Term"]
    ,s=go['%']*30,c=np.log10(go['P-Value']),marker='o',alpha=1,cmap=cm,label=go['%']
    ,edgecolors='black')
#position=fig.add_axes([0.15, 0.05, 0.7, 0.03])#位置[左,下,右,上]
#plt.colorbar(cax=position,orientation='horizontal')
#plt.legend(loc=2, bbox_to_anchor=(1.2,0.85), ncol=1, frameon=True, fontsize=12,
#handlelength=2,  borderpad = 1.8,
#handletextpad=1, title='enrichment', scatterpoints = 1)

#plt.legend()

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
plt.title(title,fontsize=15)
path=r"G:\huangjingying\data\starve\plots\\"
plt.savefig(path+title+'_go_mrnanotchange_rpfup+'+title+'.png',dpi=600, bbox_inches='tight')


#go['%']
#Out[174]: 
#0     7.1
#1     3.8
#2     2.2
#3     1.6
#4     1.1
#5     1.1
#6     1.7
#7     6.0
#8     2.1
#9     1.1
#10    1.2
#11    1.1
#####--plot legand
plt.figure(figsize=(4,8))
ax = plt.subplot(111)
test1=plt.scatter(go.loc[11,'Sublist'],go.loc[11,"Term"]
    ,s=go.loc[11,'%']*30,c='white',marker='o',alpha=1,label=go.loc[11,'%']
    ,edgecolors='black')

n=0
test2=plt.scatter(go.loc[n,'Sublist'],go.loc[n,"Term"]
    ,s=go.loc[n,'%']*30,c='white',marker='o',alpha=1,label=go.loc[n,'%']
    ,edgecolors='black')
plt.legend(loc=2, bbox_to_anchor=(1.2,0.85),borderaxespad = 0.,frameon=False,ncol=2,title='enrichment(%)',fontsize=15)
plt.savefig(path+title+'_go_mrnanotchange_rpfup+'+title+'_legand.png',dpi=600, bbox_inches='tight')
