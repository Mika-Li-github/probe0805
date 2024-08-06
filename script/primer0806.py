# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:46:26 2024

@author: mikali
"""
from Bio.SeqUtils import MeltingTemp as mt

import os
import pandas as pd

os.chdir('C:/Users/mikali/Desktop/githouse/probe0805/data')
seq381 = pd.read_csv('seq381.csv')


# 使用最近邻法计算Tm
tm_nn = mt.Tm_NN(seq381['Variant sequences'].values[0])

seq381['tm_nn'] = seq381['Variant sequences'].apply(mt.Tm_NN)

#我現在要來拆兩對了
#我覺得要補ref，minor, 還有把沒序列的去掉
#biopython似乎做得到?
#可去來源看看
#rs116046827

from Bio import Entrez
import xml.etree.ElementTree as ET
import re

#b5eca0fdc6438dd7367e13025cf005045d08	
Entrez.email = "akispen01@gmail.com"
Entrez.api_key = "b5eca0fdc6438dd7367e13025cf005045d08"


# fetch info
#rsid = 'rs12769205'  # multiple alt bases
rsid = "rs116046827"  # single alt base
handle = Entrez.efetch(db="SNP", id=rsid, retmode="text")
# result is in xml format
xml_str = handle.readline().strip()

# parse the xml string
myroot = ET.fromstring(xml_str)

# parsing from DOCSUM
# SEQ=[G/A/...]  The first character is the reference base,
# followed by alternative bases
docsum_txt = myroot.find('DOCSUM').text
print(docsum_txt)
# regex to punch out A/G/... in memory parentheses
ptn = re.compile(r'SEQ=\[(.+)\]')
bases = ptn.search(docsum_txt).groups()[0].split('/')  # assume match

# print in "ref>alt" format
print(', '.join([bases[0] + '>' + base for base in bases[1:]]))
#好感動，是ok的

