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
#--------testing
file='WT-1_RES_res.filter.txt'
sample=file[:9]
df=pd.read_csv(file,sep="\t",header=1)
df.columns

#--get ratio for each site
df['res_reads']=pd.DataFrame(df['AD:DP'].str.split(":",expand=True))[0].astype(int)
df['all_reads']=pd.DataFrame(df['AD:DP'].str.split(":",expand=True))[1].astype(int)
df['res_ratio']=df['res_reads']/df['all_reads']

#--get sum ratio of each 1Mbp
df['length']=1000000
df['ranges']=df['Start(0base)']//df['length']
df['ranges'].drop_duplicates()
