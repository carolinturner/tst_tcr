## TCRseq_Table_script: Make Table 2

# download Supplementary File S1 from paper and save as 'FileS1.csv'

library(tidyverse)
library(gtsummary)

meta <- read.csv("data/FileS1.csv")

tcrseq_table <- meta %>%
  filter(TCRseq == "Yes") %>%
  select(Sample,Sex,Age_years,Ethnicity) %>%
  rename("Age"=Age_years)

tbl_summary(tcrseq_table, by="Sample", missing = "ifany") %>% 
  as_gt() %>% 
  gt::gtsave("data/Table2_TCRseq_full.docx")