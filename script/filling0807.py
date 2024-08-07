# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 09:18:19 2024

@author: mikali
"""

#補補位點資訊
import os

# 打印当前工作目录
print("当前工作目录:", os.getcwd())

# 更改工作目录
new_directory = 'C:/Users/mikali/Desktop/githouse/probe0805/data'  # 替换为你想要的目录路径
os.chdir(new_directory)

import pandas as pd

seq381 = pd.read_csv('seq381.csv')

#沒有序列的資料可能是merge了，所以先更改其dbSNP ID，交互作用因為三個都在表中，可以先去掉，不過要記錄
#先去掉交互作用
'''
	dbSNP ID	Chromosome	Position
383	rs30187 x rs116488202	5 x 6	96788627 x 31377139
384	rs10045403 x rs116488202	5 x 6	96812030 x 31377139
'''
# 删除索引 383 和 384
seq381 = seq381.drop(index=[383, 384])
