# TCRseq_script_12: identify most public and most abundant metaclones in D7 TST and Mtb-reactive T cell data from Musvosvi et al.

library(tidyverse)

## Step 1: D7 TST data ####
d1 <- read.csv("data/metaclonotypist_results_down-sampled_beta_expanded_gr0.csv") %>%
  filter(tissue == "TST_D7") %>%
  group_by(mc.index) %>%
  mutate(publicity = length(unique(sample)),
         count_by_mc = sum(duplicate_count)) %>%
  ungroup()

# most public: metaclone 8 (look up in File S5 = cluster 57)
ind1 <- d1 %>% filter(publicity == max(publicity)) %>% pull(mc.index) %>% unique()
# most abundant: metaclone 21 (look up in File S5 = cluster 6)
ind2 <- d1 %>% filter(count_by_mc == max(count_by_mc)) %>% pull(mc.index) %>% unique()

## Step 2: data for Mtb-reactive TCRs (Musvosvi et al., Nat Med, 2023) ####
dat <- read.csv("data/Musvosvi_TableS2.csv",row.names = 1) %>%
  filter(flag == "GOOD") %>%
  select(batch_bc_well,donorId,Cell.Type,Vb,Jb,CDR3b) %>%
  na.omit() %>%
  mutate(bioidentity = paste0(Vb,CDR3b,Jb))
length(unique(dat$donorId)) # 70 participants
length(unique(dat$batch_bc_well)) # 21,212 cells

# define metaclonotypist search pattern
mc <- read.csv("data/FileS4.csv")
mc_regex <- mc %>% select(regex,Vs,index) %>% unique()

# metaclone matches
summary <- data.frame()
# loop through all regex
for (i in (1:nrow(mc_regex))){
  # define search patterns and metaclone index
  regex <- mc_regex[i,1]
  V <- mc_regex[i,2]
  index <- mc_regex[i,3]
  # look for matches in combined data
  print(paste0("checking metaclonotypist index ",i," of ",nrow(mc)))
  match <- dat %>% filter(str_detect(bioidentity,regex) & str_detect(bioidentity,V))
  # add columns to indicate match and metaclone index
  if (dim(match)[1] != 0){
    match$mc.index <- index
    match$mc.match <- 1
    # add to summary file
    summary <- rbind(summary,match)
  }
}
# check most public/abundant metaclones
results <- summary %>%
  group_by(mc.index) %>%
  mutate(publicity = length(unique(donorId)),
         abundance = length(unique(batch_bc_well))) %>%
  ungroup()
most_public <- results %>%
  arrange(desc(publicity)) %>%
  select(mc.index,publicity,abundance) %>%
  slice_head(n=1) # index 156 (found in 14 participants), look up in File S5 = cluster 6171
most_abundant <- results %>%
  arrange(desc(abundance)) %>%
  select(mc.index,publicity,abundance) %>%
  slice_head(n=1) # index 100 (matching 47 cells), look up in File S5 = cluster 7654
