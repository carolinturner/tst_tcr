# TCRseq_script_4: Calculate abundance of TCR alpha sequences associated with donor-unrestricted T cells

library(data.table)
library(tidyverse)

# load data
a.full <- fread("data/combined_alpha.csv.gz")
a.down <- fread("data/combined_subsampled_alpha.csv.gz")

# MAITs: TRAV1-2-TRAJ33/12/20
mait.full <- a.full %>%
  filter(v_call == "TRAV1-2" & (j_call %in% c("TRAJ33","TRAJ12","TRAJ20"))) %>%
  group_by(tissue,UIN) %>%
  summarise(inv.count = sum(duplicate_count)) %>%
  ungroup()
mait.down <- a.down %>%
  filter(v_call == "TRAV1-2" & (j_call %in% c("TRAJ33","TRAJ12","TRAJ20"))) %>%
  group_by(tissue,UIN) %>%
  summarise(inv.count = sum(duplicate_count)) %>%
  ungroup()

# GEMs: TRAV1-2-TRAJ9
gem.full <- a.full %>%
  filter(v_call == "TRAV1-2" & j_call == "TRAJ9") %>%
  group_by(tissue,UIN) %>%
  summarise(inv.count = sum(duplicate_count)) %>%
  ungroup()
gem.down <- a.down %>%
  filter(v_call == "TRAV1-2" & j_call == "TRAJ9") %>%
  group_by(tissue,UIN) %>%
  summarise(inv.count = sum(duplicate_count)) %>%
  ungroup()

# iNKTs: TRAV10-TRAJ18
nkt.full <- a.full %>%
  filter(v_call == "TRAV10" & j_call == "TRAJ18") %>%
  group_by(tissue,UIN) %>%
  summarise(inv.count = sum(duplicate_count)) %>%
  ungroup()
nkt.down <- a.down %>%
  filter(v_call == "TRAV10" & j_call == "TRAJ18") %>%
  group_by(tissue,UIN) %>%
  summarise(inv.count = sum(duplicate_count)) %>%
  ungroup()

# calculate percentage
dat.full <- a.full %>%
  group_by(tissue,UIN) %>%
  summarise(total.count = sum(duplicate_count)) %>%
  ungroup()
dat.down <- a.down %>%
  group_by(tissue,UIN) %>%
  summarise(total.count = sum(duplicate_count)) %>%
  ungroup()

mait.pct.full <- dat.full %>%
  left_join(mait.full) %>%
  replace(is.na(.), 0) %>%
  mutate(inv.pct = inv.count/total.count*100,
         dataset = "full",
         type = "MAIT")
mait.pct.down <- dat.down %>%
  left_join(mait.down) %>%
  replace(is.na(.), 0) %>%
  mutate(inv.pct = inv.count/total.count*100,
         dataset = "down",
         type = "MAIT")

gem.pct.full <- dat.full %>%
  left_join(gem.full) %>%
  replace(is.na(.), 0) %>%
  mutate(inv.pct = inv.count/total.count*100,
         dataset = "full",
         type = "GEM")
gem.pct.down <- dat.down %>%
  left_join(gem.down) %>%
  replace(is.na(.), 0) %>%
  mutate(inv.pct = inv.count/total.count*100,
         dataset = "down",
         type = "GEM")

nkt.pct.full <- dat.full %>%
  left_join(nkt.full) %>%
  replace(is.na(.), 0) %>%
  mutate(inv.pct = inv.count/total.count*100,
         dataset = "full",
         type = "iNKT")
nkt.pct.down <- dat.down %>%
  left_join(nkt.down) %>%
  replace(is.na(.), 0) %>%
  mutate(inv.pct = inv.count/total.count*100,
         dataset = "down",
         type = "iNKT")

summary <- rbind(mait.pct.full,mait.pct.down,gem.pct.full,gem.pct.down,nkt.pct.full,nkt.pct.down)
write.csv(summary,"data/invariant_alpha_tcrs.csv",row.names = F)