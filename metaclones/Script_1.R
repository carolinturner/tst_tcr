# Script_1 (metaclone discovery): re-formatting of HLA data

# save Table S3 from the manuscript as TableS3.csv in a sub-directory called data

library(tidyverse)

# load data
dat <- read.csv("data/TableS3.csv")

# select relevant columns and split HLA gene and allele
dat1 <- dat %>%
  select(UIN, HLA) %>%
  separate_wider_delim(cols = HLA, delim = "_", names = c("HLA","allele")) %>%
  select(-HLA) %>%
  separate_wider_delim(cols = allele, delim = "*", names = c("HLA_gene","allele2"), cols_remove = FALSE) %>%
  select(-allele2) %>%
  separate_wider_delim(cols = allele, delim = ":", names = c("allele", "field2")) %>%
  select(-field2)

# add allele duplicate where second allele not specified
d1 <- dat1 %>%
  group_by(UIN,HLA_gene) %>%
  filter(n() == 1) %>%
  slice(rep(1:n(), each =2)) %>%
  ungroup()
d2 <- dat1 %>%
  group_by(UIN,HLA_gene) %>%
  filter(n() > 1) %>%
  ungroup()
dat2 <- rbind(d1,d2)

# pivot wider
dat3 <- dat2 %>%
  group_by(HLA_gene,UIN) %>%
  mutate(HLA_group = row_number(),
         HLA_id = paste0(HLA_gene,".",HLA_group)) %>%
  ungroup() %>%
  select(-HLA_gene,-HLA_group) %>%
  pivot_wider(id_cols = "UIN", names_from = "HLA_id", values_from = "allele") %>%
  column_to_rownames("UIN")

# combine DQA and DQB alleles (all four combinations)
dat4 <- dat3 %>%
  mutate(DQAB1 = paste0(DQA1.1,"_",DQB1.1),
         DQAB2 = paste0(DQA1.1,"_",DQB1.2),
         DQAB3 = paste0(DQA1.2,"_",DQB1.1),
         DQAB4 = paste0(DQA1.2,"_",DQB1.2)) %>%
  select(-DQA1.1,-DQA1.2,-DQB1.1,-DQB1.2)

# set NA values (where no HLA allele imputation is provided) to MISSING
dat5 <- replace(dat4, is.na(dat4),"MISSING")
  
# save to file
write.csv(dat5,"data/hladata.csv",row.names = T)
