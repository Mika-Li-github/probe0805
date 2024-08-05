# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 16:17:29 2024

@author: mikali
"""

import pandas as pd

from pandasgwas.get_associations import get_associations_by_efo_id
associations = get_associations_by_efo_id('EFO_0003898')


# 按照A欄位合併表格
panda_A = pd.merge(associations.strongest_risk_alleles, associations.associations, on='associationId', how='outer')

# 使用 str.rsplit() 以最右邊的 '-' 分割 A 欄位
panda_A[['rsId', 'riskAllele']] = panda_A['riskAlleleName'].str.rsplit('-', n=1, expand=True)

import os

new_dir = "C:/Users/mikali/Desktop/githouse/probe0805/data"
os.chdir(new_dir)
print(f"Changed directory to: {os.getcwd()}")

#index=False可以避免匯出行索引。
panda_A.to_csv('panda_A.csv', index=False)

