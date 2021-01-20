#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 14:51:32 2019

@author: Yufan A. Liu
"""

# script of contrast plotting.

import pandas as pd
import sys
import seaborn as sns
import warnings
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def help_doc():
    print("")
    # 输入说明
    print("--------------------------------------------------------------------------")
    print("     Users helping documentation, support 2 files for contrast.           ")
    print("--------------------------------------------------------------------------")
    print("      Parameter                        Description                        ")
    print("         -1              Output file 1 from SPRINT software(aimed file)   ")
    print("         -2              Output file 2 from SPRINT software(contrast file)")
    print("         -g                     Gene name or Gene symbol                  ")
    print("         -id1                First sample ID tag for plotting             ")
    print("         -id2                Second sample ID tag for plotting            ")
    print("         -gtf                   GTF format file end with .gtf             ")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("Connect to liuyf19@mails.tsinghua.edu.cn when questions arise.            ")
    print("Test version. Every parameter is necessary and absolute path is required.")
    print("--------------------------------------------------------------------------")
    print("")
    sys.exit(0)


# 读取各个参数
if len(sys.argv) < 2:
    warnings.warn("ALL PARAMETERS are required! ")
    help_doc()


def main():
    i = 0
    while i < len(sys.argv):
        if sys.argv[i] == '-1':
            f1 = sys.argv[i + 1]
        elif sys.argv[i] == '-2':
            f2 = sys.argv[i + 1]
        elif sys.argv[i] == '-g':
            gene_selected = sys.argv[i + 1]
        elif sys.argv[i] == '-id1':
            id1 = sys.argv[i + 1]
        elif sys.argv[i] == '-id2':
            id2 = sys.argv[i + 1]
        elif sys.argv[i] == '-gtf':
            gtfReadFile = sys.argv[i + 1]
        i = i + 1

    '''
    #调试用
    bed_file = 'gene.fil.bed'
    f1 = 'GOU3041-3_RES_res.filter.txt'
    f2 = 'WT.txt'
    gene_selected = 'ver-1'
    id1 = 'GOU3041'
    id2 = 'WT'
    '''

    gtfFile = pd.read_table(gtfReadFile, sep="\t", skiprows=5, header=None)
    gtfFile = gtfFile[gtfFile.iloc[:, 2].str.contains('gene')]
    gtfRaw = gtfFile.iloc[:, [0, 3, 4, 6]]
    gtfRaw.columns = ['ChromType', 'Start', 'End', 'StrandType']
    gtfSplit = pd.DataFrame(gtfFile.iloc[:, [8]].iloc[:, 0].str.split(r";|\"", expand=True)).iloc[:, [1, 7, 10, 13]]
    gtfSplit.columns = ['GeneID', 'GeneName', 'DBSource', 'GeneBiotype']
    gtfBedFile = gtfRaw.join(gtfSplit).reset_index(drop=True)


    # 处理bed文件
    locus_file = gtfBedFile
    gene_name = list(locus_file.GeneName)
    start_site = locus_file.Start
    stop_site = locus_file.End

    # 处理输入文件列表，使用1-base的数据
    header = [id1, id2]
    f1_table = pd.read_table(f1, sep="\t")
    f1_table['ratio'] = round(f1_table['AD'] / f1_table['DP'], 2)
    tb1 = f1_table.drop(['Start(0base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis=1)
    f2_table = pd.read_table(f2, sep="\t")
    f2_table['ratio'] = round(f2_table['AD'] / f2_table['DP'], 2)
    tb2 = f2_table.drop(['Start(0base)', 'Supporting_reads', 'Strand', 'Type', 'AD', 'DP'], axis=1)

    # 提取gene所在的染色体位置

    pos = gene_name.index(gene_selected)
    chrom = 'chr' + str(locus_file.iloc[pos, 0])
    start = start_site[pos]
    stop = stop_site[pos]

    tb1_gene = tb1[(tb1['#Chrom'] == chrom) & (tb1['End(1base)'] >= start) & (tb1['End(1base)'] <= stop)].drop(
        ['#Chrom'], axis=1)
    tb2_gene = tb2[(tb2['#Chrom'] == chrom) & (tb2['End(1base)'] >= start) & (tb2['End(1base)'] <= stop)].drop(
        ['#Chrom'], axis=1)

    # 获取绘图列表并进行索引填充
    tb1_gene = tb1_gene.set_index(['End(1base)'])
    tb2_gene = tb2_gene.set_index(['End(1base)'])

    for i in range(start, stop):
        if i not in tb1_gene.index:
            tb1_gene.loc[i] = 0
        if i not in tb2_gene.index:
            tb2_gene.loc[i] = 0
    tb1_gene.sort_index(inplace=True)
    tb2_gene.sort_index(inplace=True)

    # 合并
    heat_table = pd.merge(tb1_gene, tb2_gene, left_index=True, right_index=True)
    heat_table.columns = header
    heat_table = heat_table.T

    # 绘图
    fig, ax = plt.subplots(figsize=(20, 2.5))
    font = {'family': 'Times New Roman', 'Weight': 'normal', 'Size': 15, }
    sns.heatmap(heat_table, cmap=colors.LinearSegmentedColormap.from_list('Mycmap',
                                                                          ['#B4B4B4', '#FF0000', '#FFA500', '#FFFF00',
                                                                           '#00FF00', '#007FFF', '#0000FF', '#8B00FF']),
                vmax=1.0, vmin=0)

    ax.set_title(gene_selected + " editing ratio graph." + chrom + " positon:" + str(start) + ":" + str(stop), font)
    ax.set_xlabel("Position", font)
    ax.set_ylabel("Sample ID", font)
    # ax.set_yticks(range(-0.5, 1.5))

    #  tick_locator = ticker.MaxNLocator(6)
    #  ax.xaxis.set_major_locator(tick_locator) #可以显示的最多的横轴坐标数目
    ax.set_yticklabels(header, va="center", rotation=45)
    # ax.set_xticklabels(rotation = 45)
    # plt.xticks(range(start, stop, int((stop-start)/5)))
    # ax.set_xticks(np.arange(start, stop, int((stop-start)/5)))
    # ax.set_xticks([])
    ax.axhline(1, linewidth=2, c='black')

    fig.tight_layout()

    fig.savefig(gene_selected + ".pdf")




main()
