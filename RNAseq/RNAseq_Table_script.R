## RNAseq_Table_script: Make Table 1

# download Supplementary File S1 from paper and save as 'FileS1.csv'

library(tidyverse)
library(gtsummary)

meta <- read.csv("data/FileS1.csv")

rnaseq_table <- meta %>%
  filter(RNAseq == "Yes") %>%
  select(Sample,Sex,Age_years,Ethnicity) %>%
  rename("Age"=Age_years)

tbl_summary(rnaseq_table, by="Sample", missing = "ifany") %>% 
  as_gt() %>% 
  gt::gtsave("data/Table1_RNAseq.docx")