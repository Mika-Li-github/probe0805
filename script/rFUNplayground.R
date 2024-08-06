new_values <- tibble(
  `dbSNP ID` = "rs116488202",
  Chromosome = 6,
  Position = 31377139,
  `ref allele` = NA,
  "minor allele (Alternative)" = NA,
  `minor allele` = NA,
  `MAF`  = NA#,
  # `Variant sequences` = NA,
  # Extracted = NA,
)


merged_table <- rbind(snp388, new_values)
merged_table <- merged_table %>% filter(Chromosome != "5 x 6")
snp_sequence_i <- getBM(attributes = c("refsnp_id", 
                                       "snp"), 
                        filters = c("snp_filter", "chr_name", "start", "upstream_flank", "downstream_flank"), 
                        values = list(merged_table$`dbSNP ID`, merged_table$Chromosome, merged_table$Position, 800, 800), 
                        mart = ensembl,
                        checkFilters = FALSE, 
                        bmHeader = TRUE
)
#為甚麼還是只有381個
