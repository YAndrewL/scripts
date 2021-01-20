# -*- coding: utf-8 -*-                                                                  
# @Time    : 2020/11/6 11:16 上午                                                          
# @Author  : Yufan A. Liu                                                                
                                                                                         
import re                                                                                

# 去掉高频的突变 AO/DP > 0.2                                                                                         
with open("Light.lowfreq.snp.vcf", 'w') as out:                                          
    with open("Light.vcf", 'r') as f:                                                    
        lines = f.readlines()                                                            
        for lin in lines:                                                                
            if lin[0] == '#':                                                            
                out.write(lin)                                                           
            else:                                                                        
                try:                                                                     
                    AO = re.compile(r'AO=.*?;').search(lin).group()                      
                    AO_value = int(AO[AO.index("=") + 1: AO.index(";")])                 
                    DP = re.compile(r'DP=.*?;').search(lin).group()                      
                    DP_value = int(DP[DP.index("=") + 1: DP.index(";")])                 
                    if AO_value / DP_value < 0.2 and "TYPE=snp" in lin:                  
                        out.write(lin)                                                   
                except:                                                                  
                    pass                                                                 
# 注释的merge                                                                                         
vcf_file = pd.read_table("Dark_genes.txt", sep="\t", header=0, skiprows=1, low_memory=False)
anno_tb = pd.read_excel("Creinhardtii_281_v5.5.mergedAllID_withAnnotation.xlsx", header=0)

x = pd.merge(vcf_file, anno_tb, left_on=['GeneId'], right_on=['ID_CRL'], how='left')
x.to_csv("Dark.fil.ann.txt", sep='\t', index = 0)