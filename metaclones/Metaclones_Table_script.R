## Metaclones_table_script: Display selected metaclones in table

library(tidyverse)

## Table 3A: most public metaclones per HLA class II allele in metaclone discovery dataset
# get lead HLA association per metaclone
dat <- read.csv("data/FileS5.csv") %>%
  group_by(index) %>%
  slice_min(pvalue,n=1) %>%
  ungroup() %>%
  mutate(HLApos = paste0(count_allele,"/",total_allele),
         HLAneg = paste0(count_other,"/",total_other)) 
# get most public metaclone per HLA allele
dat2 <- dat %>%
  arrange(desc(count_allele)) %>%
  group_by(hla) %>%
  slice_head(.,n=1) %>%
  ungroup()
# use top 10 most public as Table 3A
tab3A <- dat2 %>%
  arrange(desc(count_allele)) %>%
  slice_head(.,n=10) %>%
  select(index,consensus,Vs,hla,HLApos,HLAneg,pvalue,odds_ratio)
tab3A$odds_ratio <- round(tab3A$odds_ratio,digits =1)
write.csv(tab3A,"data/Table3A.csv",row.names = F)

## Table 3B: most public metaclones in down-sampled TST D7 dataset (n=128)
# get metaclone hits and calculate publicity
df <- read.csv("data/metaclonotypist_results_down-sampled_beta_expanded_gr0.csv") %>%
  filter(tissue == "TST_D7") %>%
  group_by(mc.index) %>%
  mutate(publicity1 = length(unique(sample))) %>%
  ungroup() %>%
  select(mc.index,publicity1) %>%
  unique() %>%
  dplyr::rename(index = "mc.index")
# annotate with metaclone description
df2 <- left_join(df,dat) %>%
  mutate(publicity = paste0(publicity1,"/",128)) %>%
  select(index,publicity1,publicity,consensus,Vs,hla,pvalue,odds_ratio)
# use top 10 most public as Table 3B
tab3B <- df2 %>%
  arrange(desc(publicity1)) %>%
  slice_head(.,n=10) %>%
  select(-publicity1)
write.csv(tab3B,"data/Table3B.csv",row.names = F)