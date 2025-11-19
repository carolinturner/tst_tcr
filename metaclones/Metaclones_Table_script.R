## Metaclones_table_script: Display selected metaclones in table

library(tidyverse)

## Table 3: most public metaclones in down-sampled TST D7 dataset (n=128)
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
tab3 <- df2 %>%
  arrange(desc(publicity1)) %>%
  slice_head(.,n=10) %>%
  select(-publicity1)
write.csv(tab3,"data/Table3.csv",row.names = F)