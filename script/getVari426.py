# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 09:27:08 2024

@author: mikali
"""


from pandasgwas.get_variants import get_variants_by_efo_id
snps = get_variants_by_efo_id('EFO_0003898')

panda_S = snps.variants


# 从字典列表中提取信息
panda_S['chromosomeName'] = panda_S['locations'].apply(lambda locs: [loc['chromosomeName'] for loc in locs])
panda_S['chromosomePosition'] = panda_S['locations'].apply(lambda locs: [loc['chromosomePosition'] for loc in locs])
panda_S['region'] = panda_S['locations'].apply(lambda locs: [loc['region'] for loc in locs])
panda_S['regionName'] = panda_S['region'].apply(lambda locs: [loc['name'] for loc in locs])


# 打印结果
print(panda_S)

#接下來要處理panda_S其他欄
# 定義一個函數來過濾符合條件的 dict
# 定義一個函數來篩選字典
def filter_dicts(dict_list):
    return [gene for gene in dict_list if gene['source'] == 'Ensembl' and gene['distance'] == 0]

# 使用 apply 方法來應用篩選函數
panda_S['filtered_genomicContexts'] = panda_S['genomicContexts'].apply(filter_dicts)

print("\n篩選後的 panda_S:")
print(panda_S)


# 使用 apply 方法來篩選列表長度大於 1 的行
more_than_1_gene_panda_S = panda_S[panda_S['filtered_genomicContexts'].apply(len) > 1]

print("\n篩選後的 panda_S:")
print(more_than_1_gene_panda_S)

# 定義一個函數來提取字符
def extract_char(value):
    # 檢查是否為空列表
    if value == "[]":
        return None  # 或者返回其他標記，例如 "N/A"
    # 去掉方括號和單引號，然後返回第一個字符
    return value.strip("[]'")

# 使用 apply 方法來提取字符
panda_S['extracted_chrName'] = panda_S['chromosomeName'].apply(extract_char)

print("\n提取後的 DataFrame:")
print(panda_S)

# 定義一個函數來提取字符
def extract_char(value):
    # 檢查列表是否為空
    if len(value) == 0:
        return None  # 或者返回其他標記，例如 "N/A"
    # 返回列表中的第一個元素
    return str(value[0])

# 使用 apply 方法來提取字符
panda_S['extracted_chrName'] = panda_S['chromosomeName'].apply(extract_char)
panda_S['extracted_chrPos'] = panda_S['chromosomePosition'].apply(extract_char)
panda_S['extracted_region'] = panda_S['regionName'].apply(extract_char)


print("\n提取後的 DataFrame:")
print(panda_S)

panda_S_dropped = panda_S.drop(columns=['locations', 'chromosomeName', 'chromosomePosition', 'region'])
panda_S_dropped.to_csv('panda_S.csv', index=False)

