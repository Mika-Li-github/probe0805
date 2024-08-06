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
# listDatasets(mart)
# searchDatasets(mart, pattern = "hsa")
# searchFilters(mart, pattern = "ref")

ensembl <- useEnsembl(biomart = "snps")
ensembl <- useDataset(dataset = "hsapiens_snp", mart = ensembl)
snp388 <- read_csv("SNP388.csv")

snp_sequence <- getBM(attributes = c("refsnp_id", 
                                     "snp"), 
                      filters = c("snp_filter", "chr_name", "start", "upstream_flank", "downstream_flank"), 
                      values = list(snp388$`dbSNP ID`, snp388$Chromosome, snp388$Position, 800, 800), 
                      mart = ensembl,
                      checkFilters = FALSE, 
                      bmHeader = TRUE
                      )



#檢查一下有沒有重複的
value_counts <- snp_sequence %>% count(`Variant name`) %>% 
  filter(n > 1)
#呵呵，沒有

# 根據B和D進行左連接
merged_table <- snp388 %>%
  left_join(snp_sequence, by = c("dbSNP ID" = "Variant name"))

# 使用正则表达式提取"%%"之间的字符
merged_table$Extracted <- str_extract(merged_table$`Variant sequences`, "%(.*?)%")

# 去掉提取结果中的"%"
merged_table$Extracted <- gsub("%", "", merged_table$Extracted)

write_csv(merged_table, "seq381.csv")
#然後就去設計探針吧

#來用chromosome找剩下序列?
# 找出B變數為NA的整列觀測值
rows_with_na_in_B <- merged_table[is.na(merged_table$`Variant sequences`), ]

#我把interaction那個沒有的也加入了 #不，事實上本來就有
# 定义要添加的值，并将其他变量填充为NA
# new_values <- tibble(
#   `dbSNP ID` = "rs116488202",
#   Chromosome = 6,
#   Position = 31377139,
#   `ref allele` = NA,
#   "minor allele (Alternative)" = NA,
#   `minor allele` = NA,
#   `MAF`  = NA,
#   `Variant sequences` = NA,
#   Extracted = NA,
# )
# 
# # 将新值添加到数据框
# merged_table <- rbind(merged_table, new_values)
# #也把紀錄interaction的觀測值刪了
# merged_table <- merged_table %>% filter(Chromosome != "5 x 6")
# #再找一次序列
# snp_sequence_i <- getBM(attributes = c("refsnp_id", 
#                                      "snp"), 
#                       filters = c("snp_filter", "chr_name", "start", "upstream_flank", "downstream_flank"), 
#                       values = list(merged_table$`dbSNP ID`, merged_table$Chromosome, merged_table$Position, 800, 800), 
#                       mart = ensembl,
#                       checkFilters = FALSE, 
#                       bmHeader = TRUE
# )
# #開補序列
# # rows_with_na_in_B <- merged_table[is.na(merged_table$`Variant sequences`), ]
# # 
# # snp_sequence_na <- getBM(attributes = c("refsnp_id", 
# #                                      "snp"), 
# #                       filters = c("chr_name", "start", "upstream_flank", "downstream_flank"), 
# #                       values = list(rows_with_na_in_B$Chromosome, rows_with_na_in_B$Position, 800, 800), 
# #                       mart = ensembl,
# #                       checkFilters = FALSE, 
# #                       bmHeader = TRUE
# # )
# #跑超久是怎樣(11:51開始，11:55都還沒結束)
# #先去python看?下一步用的新工具也許可以補這個功能
# #當前發現紀錄於mindmap「例外」中。
# #想要用序列找，若是用getSequence()，會需要end。但是又似乎沒有snp的序列資料
# #而其他還有merged的rsId，位置不同，有爭議，先不要動
# 
# #先把沒爭議的序列做成表格
# 
# snp_sequence_inte <- getBM(attributes = c("refsnp_id", 
#                                        "snp"), 
#                         filters = c("snp_filter", "chr_name", "start", "upstream_flank", "downstream_flank"), 
#                         values = list("rs116488202", 6, 31377139, 800, 800), 
#                         mart = ensembl,
#                         checkFilters = FALSE, 
#                         bmHeader = TRUE
# )
