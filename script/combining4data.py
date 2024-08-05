# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:29:47 2024

@author: mikali
"""

import os

new_dir = "C:/Users/mikali/Desktop/githouse/probe0805/data"
os.chdir(new_dir)

import pandas as pd

ensembl_export_EFO0003898 = pd.read_csv('ensembl_export_EFO0003898.csv')
panda_A = pd.read_csv('panda_A.csv')
panda_S = pd.read_csv('panda_S.csv')
gwas = pd.read_csv("gwas-association-downloaded_2024-08-02-EFO_0003898.tsv", sep="\t")


#先來刪除沒有在人類genome裡的
filtered_ensembl_export_EFO0003898 = ensembl_export_EFO0003898[ensembl_export_EFO0003898['Genomic location (strand)'].str.match(r'^[0-9]|^X')]

# 定義一個函數來處理觀測值
def process_observation(obs):
    parts = obs.split(':')
    chromosome = str(parts[0])
    
    start_end = parts[1].split('-')
    start = str(start_end[0])
    end = str(start_end[1][:-1])  # 去除最後的+或-
    
    sign = str(start_end[1][-1])  # 獲取最後的+或-
    
    # 明確指定列名
    return pd.Series([chromosome, start, end, sign], index=['ChromosomeName', 'Start', 'End', 'Strand'])

# 假設 filtered_ensembl_export_EFO0003898 是你的 DataFrame
result = filtered_ensembl_export_EFO0003898['Genomic location (strand)'].apply(process_observation)

# 使用 .loc 進行賦值
filtered_ensembl_export_EFO0003898.loc[:, ['ChromosomeName', 'Start', 'End', 'Strand']] = result
#ensembl搞定

gwas = gwas.dropna(subset=['CHR_ID', 'CHR_POS'])

# 將 'column_name' 列轉換為字串類型
gwas['SNP_ID_CURRENT'] = gwas['SNP_ID_CURRENT'].fillna(0).astype(int).astype(str)

panda_S['extracted_chrPos'] = panda_S['extracted_chrPos'].fillna(0).astype(int).astype(str)

panda_S = panda_S.dropna(subset=['extracted_chrName'])
panda_S_dropped = panda_S.drop(columns=['regionName'])

#先把重要數據存到data
panda_S_dropped.to_csv('panda_S.csv', index=False)
panda_A.to_csv('panda_A.csv', index=False)
gwas.to_csv('gwas.csv', index=False)
filtered_ensembl_export_EFO0003898.to_csv('ensembl_export_EFO0003898.csv', index=False)


#現在要從panda_S_dropped, gwas, filtered_ensembl_export_EFO0003898找出unique rsId

# 提取三個DataFrame中的'A'變量並合併
combined = pd.concat([panda_S_dropped['rsId'], gwas['SNPS'], filtered_ensembl_export_EFO0003898['Name(s)']])

# 找出唯一值
unique_values = combined.unique()

# 將唯一值建立為新的DataFrame
fourInOne = pd.DataFrame({'A': unique_values})

# 顯示新的DataFrame
print("新的DataFrame:")
print(fourInOne)

#現在來完成必要表格(ness)
#dbSNP ID	Chromosome	Position	ref allele	minor allele (Alternative)	minor allele	MAF
#這些表格都只有dbSNP ID	Chromosome	Position，但gwas的allele可以納入參考
#首先去挖panda_S_dropped
panda_S_ness = panda_S_dropped[['rsId', 'extracted_chrName', 'extracted_chrPos']]  # 提取'A'和'B'欄位
#有重複嗎
panda_S_ness = panda_S_ness.drop_duplicates()
#有重複歐

#gwas呢
gwas_ness = gwas[['SNPS', 'CHR_ID', 'CHR_POS']]  # 提取'A'和'B'欄位
gwas_ness = gwas_ness.drop_duplicates()
#有重複歐

#filtered_ensembl_export_EFO0003898呢
ensembl_ness = filtered_ensembl_export_EFO0003898[['Name(s)', 'ChromosomeName', 'Start']]  # 提取'A'和'B'欄位
ensembl_ness = ensembl_ness.drop_duplicates()

#來統一一下名稱
panda_S_ness = panda_S_ness.rename(columns={'rsId': 'dbSNP ID', 'extracted_chrName': 'Chromosome', 'extracted_chrPos': 'Position'})
gwas_ness = gwas_ness.rename(columns={'SNPS': 'dbSNP ID', 'CHR_ID': 'Chromosome', 'CHR_POS': 'Position'})
ensembl_ness = ensembl_ness.rename(columns={'Name(s)': 'dbSNP ID', 'ChromosomeName': 'Chromosome', 'Start': 'Position'})

combined_ness = pd.concat([panda_S_ness, gwas_ness, ensembl_ness], ignore_index=True)
combined_ness = combined_ness.drop_duplicates()

#輸出一下表格
new_dir = "C:/Users/mikali/Desktop/githouse/probe0805/temp"
os.chdir(new_dir)
combined_ness.to_csv('combined_ness.csv', index=False)
#去R抓其他欄位
