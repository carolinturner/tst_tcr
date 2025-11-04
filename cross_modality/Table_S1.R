## Calculation of Spearman correlation coefficients to manually populate Table S1 

# download Supplementary File S1 from paper and save as 'FileS1.csv'

library(tidyverse)
library(rstatix)
library(data.table)

# load TST induration measurements
meta <- read.csv("data/FileS1.csv") %>%
  filter(!is.na(TST_induration1)) %>%
  mutate(TST_induration2 = ifelse(is.na(TST_induration2), 0, TST_induration2)) %>%
  rowwise() %>%
  mutate(induration = max(c_across(TST_induration1:TST_induration2))) %>%
  ungroup() %>%
  select(UIN,induration) %>%
  unique()

#### Correlation with transcription modules ####
rnaseq_meta <- read.csv("data/RNAseq_metadata.csv")
mod <- read.csv("data/Modules_Z-scores.csv") %>%
  left_join(rnaseq_meta) %>%
  select(UIN,Stimulant,module_name,module_Zscore)

# TST_D2
dat <- mod %>%
  filter(Stimulant == "TST_D2") %>%
  pivot_wider(names_from = module_name, values_from = module_Zscore)
corrdat <- dat %>% left_join(meta) %>% na.omit()
cor_test(corrdat, induration, T_cell, method = "spearman")
cor_test(corrdat, induration, CD4_T, method = "spearman")
cor_test(corrdat, induration, CD8_T, method = "spearman")
cor_test(corrdat, induration, NK, method = "spearman")
cor_test(corrdat, induration, myeloid, method = "spearman")
cor_test(corrdat, induration, CCND1, method = "spearman")

# TST_D7
dat <- mod %>%
  filter(Stimulant == "TST_D7") %>%
  pivot_wider(names_from = module_name, values_from = module_Zscore)
corrdat <- dat %>% left_join(meta) %>% na.omit()
cor_test(corrdat, induration, T_cell, method = "spearman")
cor_test(corrdat, induration, CD4_T, method = "spearman")
cor_test(corrdat, induration, CD8_T, method = "spearman")
cor_test(corrdat, induration, NK, method = "spearman")
cor_test(corrdat, induration, myeloid, method = "spearman")
cor_test(corrdat, induration, CCND1, method = "spearman")

#### Correlation with percentage of metaclones in down-sampled repertoires ####
# TST_D2
mc.pct <- read.csv("data/summary_metaclone-abundance_down-sampled_beta_expanded_gr0.csv") %>%
  filter(tissue == "TST_D2") %>%
  mutate(UIN = str_split_i(sample,"_",1)) %>%
  select(UIN,mc.pct)
corrdat <- mc.pct %>% left_join(meta) %>% na.omit()
cor_test(corrdat, induration, mc.pct, method = "spearman")

# TST_D7
mc.pct <- read.csv("data/summary_metaclone-abundance_down-sampled_beta_expanded_gr0.csv") %>%
  filter(tissue == "TST_D7") %>%
  mutate(UIN = str_split_i(sample,"_",1)) %>%
  select(UIN,mc.pct)
corrdat <- mc.pct %>% left_join(meta) %>% na.omit()
cor_test(corrdat, induration, mc.pct, method = "spearman")

# Blood
mc.pct <- read.csv("data/summary_metaclone-abundance_down-sampled_beta_expanded_gr0.csv") %>%
  filter(tissue == "Blood") %>%
  mutate(UIN = str_split_i(sample,"_",1)) %>%
  select(UIN,mc.pct)
corrdat <- mc.pct %>% left_join(meta) %>% na.omit()
cor_test(corrdat, induration, mc.pct, method = "spearman")

#### Correlation with percentage of expanded beta-chain TCRs in down-sampled repertoires ####
## Step 1: calculate percentage of expanded beta-chain TCRs
data.b <- fread("data/combined_subsampled_beta.csv.gz") 

# calculate size of repertoire per sample
summary <- data.b %>%
  group_by(tissue,sample) %>%
  summarise(total.count = sum(duplicate_count)) %>%
  ungroup()

# loop through different expansion thresholds
for (f in c(1:4)) {
  print(paste0("Expansion threshold: ",f))
  # select TCRs and calculate size of expanded repertoire per sample
  dat <- data.b %>%
    filter(duplicate_count >f) %>%
    group_by(tissue,sample) %>%
    summarise(!!paste0("count_gr_", f) := sum(duplicate_count), .groups = "drop")
  # add to summary and calculate percentage expanded TCRs
  summary <- left_join(summary,dat)
  summary <- summary %>%
    mutate(
      !!paste0("pct_gr_", f) := (!!sym(paste0("count_gr_", f)) / total.count) * 100,
      !!paste0("pct_gr_", f) := replace_na(!!sym(paste0("pct_gr_", f)), 0)
    )
}
write.csv(summary,"data/summary_pct_expanded_tcrs_down-sampled_all-tissue.csv")

## Step 2: correlation with induration
# TST_D2
exp <- read.csv("data/summary_pct_expanded_tcrs_down-sampled_all-tissue.csv", row.names = 1) %>%
  filter(tissue == "TST_D2") %>%
  mutate(UIN = str_split_i(sample,"_",1)) %>%
  select(UIN,starts_with("pct"))
corrdat <- exp %>% left_join(meta) %>% na.omit()
cor_test(corrdat, induration, pct_gr_1, method = "spearman")
cor_test(corrdat, induration, pct_gr_2, method = "spearman")
cor_test(corrdat, induration, pct_gr_3, method = "spearman")
cor_test(corrdat, induration, pct_gr_4, method = "spearman")

# TST_D7
exp <- read.csv("data/summary_pct_expanded_tcrs_down-sampled_all-tissue.csv", row.names = 1) %>%
  filter(tissue == "TST_D7") %>%
  mutate(UIN = str_split_i(sample,"_",1)) %>%
  select(UIN,starts_with("pct"))
corrdat <- exp %>% left_join(meta) %>% na.omit()
cor_test(corrdat, induration, pct_gr_1, method = "spearman")
cor_test(corrdat, induration, pct_gr_2, method = "spearman")
cor_test(corrdat, induration, pct_gr_3, method = "spearman")
cor_test(corrdat, induration, pct_gr_4, method = "spearman")

# Blood
exp <- read.csv("data/summary_pct_expanded_tcrs_down-sampled_all-tissue.csv", row.names = 1) %>%
  filter(tissue == "Blood") %>%
  mutate(UIN = str_split_i(sample,"_",1)) %>%
  select(UIN,starts_with("pct"))
corrdat <- exp %>% left_join(meta) %>% na.omit()
cor_test(corrdat, induration, pct_gr_1, method = "spearman")
cor_test(corrdat, induration, pct_gr_2, method = "spearman")
cor_test(corrdat, induration, pct_gr_3, method = "spearman")
cor_test(corrdat, induration, pct_gr_4, method = "spearman")
