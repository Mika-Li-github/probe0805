# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 14:37:11 2024

@author: mikali
"""

import os

os.chdir('C:/Users/mikali/Desktop/gitfamily/TWB/TWBv2.0_SNPs位點相關資訊')

import pandas as pd

TWBSNP位點 = pd.read_csv('TWBv2.0 SNPs位點相關資訊.csv')

TWBSNP位點.rename(columns={'dbSNP RS ID': 'rsId'}, inplace=True)

# 選取需要的變量
print(TWBSNP位點.columns)
TWBSNP = TWBSNP位點[['rsId', 'Chromosome', 'Physical Position', 'Ref Allele', 'Alt Allele', "Minor Allele", 'Minor Allele Frequency']]

os.chdir('C:/Users/mikali/Desktop/githouse/probe0805/temp')
#四大天王來的388SNP:combined_ness.csv

PandaAndWebSNP = pd.read_csv('combined_ness.csv')
print(PandaAndWebSNP.columns)

filtered_TWBSNP = TWBSNP[TWBSNP['rsId'].isin(PandaAndWebSNP['dbSNP ID'])]
#236個

filtered_TWBSNP = filtered_TWBSNP.drop_duplicates()
#事實上只有203個在四大天王裡

filtered_PandaAndWebSNP = PandaAndWebSNP[~PandaAndWebSNP['dbSNP ID'].isin(TWBSNP['rsId'])]
#四大天王中，不在TWBv2.0 SNPs位點相關資訊的，185個

os.chdir('C:/Users/mikali/Desktop/githouse/probe0805/data')
filtered_TWBSNP.to_csv('filtered_TWBSNP.csv', index=False)

#先以上面的來找序列


'''
現在要來看差補的資訊
'''
# 讀取TSV文件
tsv_file_path = 'C:/Users/mikali/Desktop/gitfamily/TWB/TWBv2.0_差補位點資訊.v2/TWB2.hg38.infoscore.v2.txt'  # 替換為你的TSV文件路徑
TWB2_hg38_infoscore_v2 = pd.read_csv(tsv_file_path, sep='\t')

'''
【不要!】存儲為CSV文件
csv_file_path = 'C:/Users/mikali/Desktop/gitfamily/TWB/TWBv2.0_差補位點資訊.v2/TWB2.hg38.infoscore.v2.csv'  # 替換為你想要的CSV文件路徑
df.to_csv(csv_file_path, index=False)

因為資訊會遺失
'''
print(TWB2_hg38_infoscore_v2.columns)

TWB2_hg38_infoscore_v2.head()

filtered_TWB2_hg38_infoscore_v2 = TWB2_hg38_infoscore_v2[TWB2_hg38_infoscore_v2['[3]ID'].isin(PandaAndWebSNP['dbSNP ID'])]
#保險起見去掉重複值
filtered_TWB2_hg38_infoscore_v2 = filtered_TWB2_hg38_infoscore_v2.drop_duplicates()
#呵呵，沒有。總共311個

#那有沒有不在那203裡面的?
filtered_TWB2_hg38_infoscore_v2_2 = filtered_TWB2_hg38_infoscore_v2[~filtered_TWB2_hg38_infoscore_v2['[3]ID'].isin(filtered_TWBSNP['rsId'])]
#歐齁，有133個
diff_from_TWB2_133 = filtered_TWB2_hg38_infoscore_v2_2[['[3]ID', '# [1]CHROM', '[2]POS', '[4]REF', '[5]ALT']]

#保險起見去掉重複值
diff_from_TWB2_133 = diff_from_TWB2_133.drop_duplicates()
#呵呵，沒有

diff_from_TWB2_133.to_csv('diff_from_TWB2_133.csv', index=False)

#388中還有哪些不在TWB
filtered_PandaAndWebSNP = filtered_PandaAndWebSNP[~filtered_PandaAndWebSNP['dbSNP ID'].isin(diff_from_TWB2_133['[3]ID'])]
#52個，酷喔
filtered_PandaAndWebSNP.to_csv('TWB2missing_SNPS.csv', index=False)

#來合併阿
#首先變數要一致
#filtered_PandaAndWebSNP
#diff_from_TWB2_133
#filtered_TWBSNP
diff_from_TWB2_133.rename(columns={'[3]ID': 'dbSNP ID', '# [1]CHROM': 'Chromosome', '[2]POS': 'Position', '[4]REF': 'ref allele', '[5]ALT': 'minor allele (Alternative)'}, inplace=True)
print(diff_from_TWB2_133.columns)

# 移除特定字元（例如 @, #, $, %）
diff_from_TWB2_133['Chromosome'] = diff_from_TWB2_133['Chromosome'].str.replace(r'[chr]', '', regex=True)

filtered_TWBSNP.rename(columns={'rsId': 'dbSNP ID', 'Physical Position': 'Position', 'Ref Allele': 'ref allele', 'Alt Allele': 'minor allele (Alternative)', 'Minor Allele': 'minor allele', "Minor Allele Frequency": 'MAF'}, inplace=True)

merged_df = pd.concat([filtered_TWBSNP, diff_from_TWB2_133, filtered_PandaAndWebSNP], ignore_index=True)

merged_df.to_csv('SNP388.csv', index=False)
