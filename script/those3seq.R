library(tidyverse)
library(biomaRt)
library(BiocManager)
setwd("C:/Users/mikali/Desktop/githouse/probe0805/data")

# ENSEMBL_MART_SNP <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# seq388 <- getSequence(chromosome = snp388$Chromosome, 
#                       start = snp388$Position,
#                       #id = snp388$`dbSNP ID`, 
#                       type = 'chromosome_name', 
#                       seqType = 'cdna', 
#                       upstream = 800, downstream = 800, mart = mart)
#  You must specify both a start and end position.
#所以我不用getSequence了
#不需要轉成FASTA，就不轉

# ?useMart
mart = useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
listDatasets(mart)
searchDatasets(mart, pattern = "hsa")
searchFilters(mart, pattern = "ref")

ensembl <- useEnsembl(biomart = "snps")
ensembl <- useDataset(dataset = "hsapiens_snp", mart = ensembl)
#snp388 <- read_csv("SNP388.csv")
snp3 <- read_csv("only3_TWBSNP.csv")

snp_sequence <- getBM(attributes = c("refsnp_id", 
                                     "snp"), 
                      filters = c("snp_filter", "chr_name", "start", "upstream_flank", "downstream_flank"), 
                      values = list(snp3$rsId, snp3$Chromosome, snp3$`Physical Position`, 800, 800), 
                      mart = ensembl,
                      checkFilters = FALSE, 
                      bmHeader = TRUE
)
# 根據B和D進行左連接
merged_table <- snp3 %>%
  left_join(snp_sequence, by = c("rsId" = "Variant name"))

# 使用正则表达式提取"%%"之间的字符
merged_table$Extracted <- str_extract(merged_table$`Variant sequences`, "%(.*?)%")

# 去掉提取结果中的"%"
merged_table$Extracted <- gsub("%", "", merged_table$Extracted)

######那其他merged的呢?#######
other3 <- tibble(
  rsId = c("rs5837881", "rs28533662", "rs72867732"),
  Chromosome = c(2), 
  `Physical Position` = c(203843041) 
)

snp_sequence <- getBM(attributes = c("refsnp_id", 
                                     "snp"), 
                      filters = c("snp_filter", "chr_name", "start", "upstream_flank", "downstream_flank"), 
                      values = list(other3$rsId, snp3$Chromosome, snp3$`Physical Position`, 800, 800), 
                      mart = ensembl,
                      checkFilters = FALSE, 
                      bmHeader = TRUE
)
# 根據B和D進行左連接
merged_table <- snp3 %>%
  left_join(snp_sequence, by = c("rsId" = "Variant name"))

# 使用正则表达式提取"%%"之间的字符
merged_table$Extracted <- str_extract(merged_table$`Variant sequences`, "%(.*?)%")

# 去掉提取结果中的"%"
merged_table$Extracted <- gsub("%", "", merged_table$Extracted)
