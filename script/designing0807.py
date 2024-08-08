# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 10:14:50 2024

@author: mikali
"""
#先往下做
#來吧，取出完整資料的
import os

# 打印当前工作目录
print("当前工作目录:", os.getcwd())

# 更改工作目录
new_directory = 'C:/Users/mikali/Desktop/githouse/probe0805/data'  # 替换为你想要的目录路径
os.chdir(new_directory)

import pandas as pd

seq381 = pd.read_csv('seq381.csv')
seq203 = seq381.dropna()
'''
	dbSNP ID	Chromosome	Position	ref allele	minor allele (Alternative)	minor allele	MAF	Extracted	Variant sequences
26	rs71559680	6	21430497	TAG	CAT	---	---	TAG/CAT	TTCTTTTTTTTTTTTTTTTTCTGTGACAGAGTCTCACTTTGCTGCCCAGGCTAGAGGGCAGTGGCACGATCTCAGCTCACTGCAAGCTTTGCCTCCTGGGTTCACGCCATTCTCCTCCCTCAGCCTCCCAAGTAGCTGGGACTACAGGTGTCCGCCACCACGCCCGGCTAATTTTTTTGTATTTTTAGTAGAGACAGGGTTTCACCGTGTTAGCCAGGATGGACTCGATCTGCTGACCTCGTGATCCGCCCATCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGAGCCTGGCCCTCTTTTTTTTTTTTTTTGAGACGGAGTTTTGCCCTTGTTGCCCAGGCTGGTGTGCAATGGTGCGATCTCAGCTCAGGGCAACCTCTGCCTTCCAGGTTCAAGAGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGATTACAGGCACGTGTCACCACACCTGGCTAATTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTGGTCAGGCTGGTCTCGAATTCCTGACCTCAGGTGATTCACCCACCTCAACCTCCCAAAATGCTGGGATTACAGGCGTGAGTCACTGTGCCCGGTCCCACCTGTATTTTCAAAATAATATTTGGCCAAATTTTTATTCCTTGAATGGAGTTATCAGGAAAAAGTAAAACAATTAATGAAATAGCATGCTTTGTGGGGAGGGAGAAAGATAATCACATGTTAATTTAACTAGGAGGTGATAAAATCCACATCTTCTATGCATAATTTATCTGACCTACTACCGACAGCATAAC%TAG/CAT%TATATAAAATATCAAAGCTTTATTTGATGCTTTCCTGGGAATTCACTATTTTACCTTCGTTGAAACTGATTAGAATATAGTACCTTTTACTGTCCTTTTAGAACATTTAGAATGAGAAGTGAGCTTGCTTTCCTAGTTTCAATTCAGAAAACTACAAGAAGTGTTTATGTAGCAGCAGCTCAGTATACATGACCCTGTGCTTAGAGATATTGGGATTAAAGCTTGGACTTCCATGATTCTTGGAGTCCATTTTAGGCATCATGTAATTCATACAAATGCCCTGTGGAGCACCTACTCGGGACCAGGTAAAGCAGGCTATTTTTACATTAGTATTTTGAGATTGCCACCTGGTAACTGTCACCTTTACTTCCTTCAGAGAAAATTTAAAAAACAAACAAACAACAACAACAAAAAAAACGGAACTCTTGTCACTGGTGTTTAATGAAGACTAATAGAAACTGATTTTCATGCAGAACAAAGAACCTAATTGGTTCCACAGTTACACAGCACTCCTAAGAAACAAACCATAGGGACTTCTAGCTAAACATGGCAAATTAAAAATATCTGTATCTATATCTATATGTATATTTCCAATCCTTCTCAAAGTCCCCCTAAGATGACAGTCAGTGGTATAGTTTGGATATTTGTCCCCCCGAATCTCATGTTAAAATGCAATCCCCAATTCTGGAGGTGAGGCAGGCAGATCACTTGAAGTTAGGAGTTCGAGACCAGGTTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAACATGATCCAGTTGGGTC
資料不完整，保險起見先不使用
'''
seq202 = seq203.drop(index = 26)
# 定义互换规则
complement = str.maketrans('ATCG', 'TAGC')

# 提取fwd和Rev
def extract_fwd_rev(sequence):
    # 分割字符串
    parts = sequence.split('%')
    fwd = parts[0]  # 第一個%左邊的部分
    rev = parts[2]  # 第二個%右邊的部分
    rev = rev.translate(complement)[::-1]  # 互換並反轉
    return fwd, rev

# 應用函數並創建新列
seq202[['fwd', 'rev']] = seq202['Variant sequences'].apply(extract_fwd_rev).apply(pd.Series)
# 將 A 和 B 列連接起來，形成新的 C 列，沒有空格
seq202['fwd_r'] = seq202['fwd'] + seq202['ref allele']
seq202['fwd_a'] = seq202['fwd'] + seq202['minor allele (Alternative)']

# 定義一個替換函數
def replace_char(char):
    replacement_dict = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    return replacement_dict.get(char, char)

# 使用 apply() 方法應用替換函數
seq202['Reversed_r'] = seq202['ref allele'].apply(replace_char)
seq202['Reversed_a'] = seq202['minor allele (Alternative)'].apply(replace_char)

seq202['rev_r'] = seq202['rev'] + seq202['Reversed_r']
seq202['rev_a'] = seq202['rev'] + seq202['Reversed_a']

#先來簡單計算一下Tm值

from Bio.SeqUtils import MeltingTemp as mt

# 示例序列
sequence = 'ACGTTGCAATGCCGTA'

# 使用Wallace法计算Tm
tm_wallace = mt.Tm_Wallace(sequence)
print(f'Tm (Wallace): {tm_wallace:.2f} °C')

# 使用GC含量法计算Tm 
tm_gc = mt.Tm_GC(sequence)
print(f'Tm (GC): {tm_gc:.2f} °C')

# 使用最近邻法计算Tm
tm_nn = mt.Tm_NN(sequence)
print(f'Tm (NN): {tm_nn:.2f} °C')

#我要採用最鄰近法

# 定義函數來計算符合條件的序列
def calculate_tm(sequences):
    results = []
    for seq in sequences:
        best_sequence = None
        best_tm = None
        closest_tm_diff = float('inf')  # 初始化為無窮大

        # 確保序列長度大於等於16
        if len(seq) >= 16:
            # 從右邊取出長度 >= 16 的序列
            for length in range(16, len(seq) + 1):
                sub_seq = seq[-length:]
                tm_value = mt.Tm_NN(sub_seq)

                # 檢查Tm值是否在範圍內
                if 61 < tm_value < 63.9:
                    # 計算與62.5的差距
                    tm_diff = abs(tm_value - 62.5)

                    # 找到最接近62.5的Tm值
                    if tm_diff < closest_tm_diff:
                        closest_tm_diff = tm_diff
                        best_sequence = sub_seq
                        best_tm = tm_value

        results.append((best_sequence, best_tm))
    return results

# 指定要計算的變數
variables_to_calculate = ['fwd_r', 'fwd_a', 'rev_r', 'rev_a']  # 你可以根據需要指定要計算的變數

# 對指定的變數進行計算
for var in variables_to_calculate:
    seq202[f'{var}_fwd_r'], seq202[f'{var}_Tm_fr'] = zip(*calculate_tm(seq202[var]))
#算了五分鐘終於算完了
# 獲取列名
column_names = seq202.columns

# 重新命名列
new_columns = {
    'fwd_r_fwd_r':'Fwd_r',
    'fwd_r_Tm_fr':'Tm_fr',
    'fwd_a_fwd_r':'Fwd_a',
    'fwd_a_Tm_fr':'Tm_fa',
    'rev_r_fwd_r':'Rev_r',
    'rev_r_Tm_fr':'Tm_rr',
    'rev_a_fwd_r':'Rev_a',
    'rev_a_Tm_fr':'Tm_ra'
    }

'''
存檔點:之後都改在seq202_Tm上面
'''
seq202_Tm = seq202.rename(columns=new_columns)


# 定義一個函數來計算字符數量並進行比較
def compare_lengths(row):
    length_A_B = len(row['Fwd_r']) + len(row['Rev_r'])  # 計算A和B的字符總長度
    length_C_D = len(row['Fwd_a']) + len(row['Rev_a'])  # 計算C和D的字符總長度
    
    if length_A_B > length_C_D:
        return pd.Series({'more': 'Fwd_r', 'probe_length': length_A_B - 1})
    else:
        return pd.Series({'more': 'Fwd_a', 'probe_length': length_C_D - 1})

# 應用函數並添加新變數到DataFrame
seq202_Tm[['more', 'probe_length']] = seq202_Tm.apply(compare_lengths, axis=1)

# 根據變數E的值選擇變數B或D的長度來計算變數C
seq202_Tm['after_allele_length'] = seq202_Tm.apply(lambda row: (len(row['Fwd_r']) if row['more'] == 'Fwd_r' else row['probe_length'] - len(row['Fwd_a'])), axis=1)


# 提取第二个%后面的字符作为变量D
def extract_d(row):
    parts = row['Variant sequences'].split('%')
    if len(parts) > 2:
        return parts[2]  # 返回第二个%后面的字符
    return ''  # 如果没有第二个%则返回空字符串

# 应用提取函数
seq202_Tm['rev_after_%2'] = seq202_Tm.apply(extract_d, axis=1)

# 根据变量E的值和变量A的值截取相应的字符串，并添加B或K的前缀
def custom_slice_with_prefix(row):
    if row['more'] == 'Fwd_a':
        substring = row['rev_after_%2'][:row['after_allele_length']]
        prefix = row['Fwd_a']  # 使用B作为前缀
    else:
        substring = row['rev_after_%2'][:row['after_allele_length']]
        prefix = row['Fwd_r']  # 使用K作为前缀
    return prefix + substring

seq202_Tm['Probe'] = seq202_Tm.apply(custom_slice_with_prefix, axis=1)

#要不，先做primer表?
'''
rs_id	primer_id	Tm	GC_content(%)	genotype_label	primer_len	SNP_chr	SNP_position	primer_sequence
'''
seq202_Tm.columns
#dbSNP ID : rs_id

#rs_id + Fwd_/Rev_ + Ref/Alt
#primer_id

#'Tm_fa/ra/fr/rr'
#Tm 

# Count(G + C)/Count(A + T + G + C) * 100%。
#GC_content(%)

#'ref allele', 'minor allele (Alternative)','Reversed_r', 'Reversed_a'
#genotype_label

#len()'Fwd_r', 'Fwd_a', 'Rev_r', 'Rev_a'
#primer_len

#Chromosome
#SNP_chr

#Position
#SNP_position

#'Fwd_r', 'Fwd_a', 'Rev_r', 'Rev_a'
#primer_sequence
primer_FR = pd.DataFrame({
    'rs_id': seq202_Tm['dbSNP ID'],
    'primer_id': seq202_Tm['dbSNP ID'].astype(str) + '_Fwd' + seq202_Tm['ref allele'],
    'Tm': seq202_Tm['Tm_fr'],
    'genotype_label': seq202_Tm['ref allele'],
    'primer_len': seq202_Tm['Fwd_r'].apply(len),  # 使用 apply 計算每行的長度
    'SNP_chr': seq202_Tm['Chromosome'],
    'SNP_position': seq202_Tm['Position'],
    'GC_content(%)': (seq202_Tm['Fwd_r'].apply(lambda x: x.count('G') + x.count('C')) / seq202_Tm['Fwd_r'].apply(len)) * 100,
    'primer_sequence': seq202_Tm['Fwd_r']
    })


primer_FA = pd.DataFrame({
    'rs_id': seq202_Tm['dbSNP ID'],
    'primer_id': seq202_Tm['dbSNP ID'].astype(str) + '_Fwd' + seq202_Tm['minor allele (Alternative)'],
    'Tm': seq202_Tm['Tm_fa'],
    'genotype_label': seq202_Tm['minor allele (Alternative)'],
    'primer_len': seq202_Tm['Fwd_a'].apply(len),  # 使用 apply 計算每行的長度
    'SNP_chr': seq202_Tm['Chromosome'],
    'SNP_position': seq202_Tm['Position'],
    'GC_content(%)': (seq202_Tm['Fwd_a'].apply(lambda x: x.count('G') + x.count('C')) / seq202_Tm['Fwd_a'].apply(len)) * 100,
    'primer_sequence': seq202_Tm['Fwd_a']
    })

primer_RR = pd.DataFrame({
    'rs_id': seq202_Tm['dbSNP ID'],
    'primer_id': seq202_Tm['dbSNP ID'].astype(str) + '_Rev' + seq202_Tm['Reversed_r'],
    'Tm': seq202_Tm['Tm_rr'],
    'genotype_label': seq202_Tm['Reversed_r'],
    'primer_len': seq202_Tm['Rev_r'].apply(len),  # 使用 apply 計算每行的長度
    'SNP_chr': seq202_Tm['Chromosome'],
    'SNP_position': seq202_Tm['Position'],
    'GC_content(%)': (seq202_Tm['Rev_r'].apply(lambda x: x.count('G') + x.count('C')) / seq202_Tm['Rev_r'].apply(len)) * 100,
    'primer_sequence': seq202_Tm['Rev_r']
    })

primer_RA = pd.DataFrame({
    'rs_id': seq202_Tm['dbSNP ID'],
    'primer_id': seq202_Tm['dbSNP ID'].astype(str) + '_Rev' + seq202_Tm['Reversed_a'],
    'Tm': seq202_Tm['Tm_ra'],
    'genotype_label': seq202_Tm['Reversed_a'],
    'primer_len': seq202_Tm['Rev_a'].apply(len),  # 使用 apply 計算每行的長度
    'SNP_chr': seq202_Tm['Chromosome'],
    'SNP_position': seq202_Tm['Position'],
    'GC_content(%)': (seq202_Tm['Rev_a'].apply(lambda x: x.count('G') + x.count('C')) / seq202_Tm['Rev_a'].apply(len)) * 100,
    'primer_sequence': seq202_Tm['Rev_a']
    })


# 合併成一個長格式的 DataFrame
primer = pd.concat([primer_FR, primer_FA, primer_RR, primer_RA], ignore_index=True)

#剩probe了:D加油
probe = pd.DataFrame({
    'rs_id': seq202_Tm['dbSNP ID'],
    'Tm': seq202_Tm['Probe'].apply(lambda x: mt.Tm_NN(x)),  # 計算每行的 Tm 值
    'probe_sequence': seq202_Tm['Probe'],
    'probe_len': seq202_Tm['Probe'].apply(len),  # 使用 apply 計算每行的長度
    'SNP_chr': seq202_Tm['Chromosome'],
    'SNP_position': seq202_Tm['Position'],
    'GC_content(%)': (seq202_Tm['Probe'].apply(lambda x: x.count('G') + x.count('C')) / seq202_Tm['Probe'].apply(len)) * 100
    })


primer.to_csv('primer.csv', index=False, encoding='utf-8')

probe.to_csv('probe.csv', index=False, encoding='utf-8')
