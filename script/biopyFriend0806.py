# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:34:43 2024

@author: mikali
"""

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

#我能否用bio python找特定RNA seq 的chromosome name, position, ref allele， minor allele, MAF，以及上下游各800 nucleotides 的序列?
#https://stackoverflow.com/questions/71912494/using-entrez-efetch-to-determine-snp-reference-allele
#有機會
from Bio import Entrez
import xml.etree.ElementTree as ET
import re

#b5eca0fdc6438dd7367e13025cf005045d08	
Entrez.email = "akispen01@gmail.com"
Entrez.api_key = "b5eca0fdc6438dd7367e13025cf005045d08"


# fetch info
#rsid = 'rs12769205'  # multiple alt bases
rsid = "rs114525117"  # single alt base
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
