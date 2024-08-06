# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 13:39:34 2024

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


#rs5837881, rs11363316, rs2066847, rs28533662, rs72867732 都不在288裡

data = {
        'rsId':['rs5837881', 'rs11363316', 'rs2066847', 'rs28533662', 'rs72867732', 'rs116488202']
        }

small_TWBSNP = TWBSNP[TWBSNP['rsId'].isin(data['rsId'])]
#只找到三個
os.chdir('C:/Users/mikali/Desktop/githouse/probe0805/data')
small_TWBSNP.to_csv('only3_TWBSNP.csv', index=False)
