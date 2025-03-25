# TCRseq_script_10: identify most public and most abundant metaclones in D7 TST

library(tidyverse)

# load data
d1 <- read.csv("data/metaclonotypist_results_down-sampled_beta_expanded_gr0.csv") %>%
  filter(tissue == "TST_D7") %>%
  group_by(mc.index) %>%
  mutate(publicity = length(unique(sample)),
         count_by_mc = sum(duplicate_count)) %>%
  ungroup()

# most public: metaclone 8 (look up in Table S4 = cluster 57)
ind1 <- d1 %>% filter(publicity == max(publicity)) %>% pull(mc.index) %>% unique()
# most abundant: metaclone 21 (look up in Table S4 = cluster 6)
ind2 <- d1 %>% filter(count_by_mc == max(count_by_mc)) %>% pull(mc.index) %>% unique()