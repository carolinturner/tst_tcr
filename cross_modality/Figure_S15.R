## Comparison of TST metrics between participants with presumed recent vs. remote exposure to TB 

# download Supplementary File S1 from paper and save as 'FileS1.csv'

library(tidyverse)
library(rstatix)
library(ggpubr)

#My_Theme
t = 10 #size of text
m = 10 #size of margin around text
tc = "black" #colour of text
My_Theme = theme(
  axis.title.x = element_text(size = t, face = "bold", margin = margin(t = m)),
  axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = 0, hjust = 0.5),
  axis.title.y.left = element_text(size = t, face = "bold", margin = margin(r = m)),
  axis.title.y.right = element_text(size = t, face = "bold", margin = margin(l = m)),
  axis.text.y = element_text(size = t, face = "bold", colour = tc),
  legend.title = element_text(size=t, face = "bold", colour = tc),
  legend.text = element_text(size=t, face = "bold", colour = tc),
  plot.title = element_text(size=t, face = "bold", colour = tc),
  strip.text = element_text(size=t, face = "bold", colour = tc),
  strip.background = element_rect(fill = "gray90", colour = "black", linewidth = 0.5),
  panel.border = element_rect(fill = NA, linewidth = 0.5, colour = tc),
  panel.background = element_rect(fill = "gray97"),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  legend.position = "right", legend.justification = "top"
)

# load indication for QUANTIferon
meta <- read.csv("data/FileS1.csv") %>%
  select(UIN, Quantiferon_indication) %>%
  filter(Quantiferon_indication %in% c("RecentContactHousehold","OccHealthScreening")) %>%
  mutate(TB.exposure = recode(Quantiferon_indication,
                              "OccHealthScreening" = "OH",
                              "RecentContactHousehold" = "HC")) %>%
  unique()

## TST induration #####
tst.ind <- read.csv("data/FileS1.csv") %>%
  filter(!is.na(TST_induration1)) %>%
  mutate(TST_induration2 = ifelse(is.na(TST_induration2), 0, TST_induration2)) %>%
  rowwise() %>%
  mutate(induration = max(c_across(TST_induration1:TST_induration2))) %>%
  ungroup() %>%
  select(UIN,induration) %>%
  unique()

dat <- inner_join(meta,tst.ind)

stats <- dat %>%
  wilcox_test(data = ., induration~TB.exposure) %>%
  add_xy_position() # n=27 vs. 86

pA <- ggplot(dat,aes(x=TB.exposure,y=induration)) +
  geom_boxplot() +
  My_Theme +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  stat_pvalue_manual(stats) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("") +
  ylab("Day 2 TST induration (mm)")
pA

## transcriptional modules ####
rnaseq_meta <- read.csv("data/RNAseq_metadata.csv")
mod <- read.csv("data/Modules_Z-scores.csv") %>%
  left_join(rnaseq_meta) %>%
  select(UIN,Stimulant,module_name,module_Zscore) %>%
  mutate(module_name = recode(
    module_name, CCND1 = "Prolif"
  ))
mod$module_name <- factor(mod$module_name,levels = c("Prolif","T_cell","CD4_T","CD8_T","myeloid","NK"))

dat <- inner_join(meta,mod) %>%
  select(UIN,TB.exposure,Stimulant,module_name,module_Zscore) 

stats <- dat %>%
  group_by(Stimulant,module_name) %>%
  wilcox_test(data = ., module_Zscore~TB.exposure) %>%
  add_xy_position() # n=30 vs. 100 on D2; n=27 vs. 67 on D7

pB <- ggplot(dat, aes(x=TB.exposure,y=module_Zscore)) +
  geom_boxplot() +
  facet_grid(Stimulant~module_name) +
  stat_pvalue_manual(stats) +
  My_Theme +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  xlab("") +
  ylab("Module Z score")
pB

## TCRseq diversity (TST_D7) ####
meta_tcr <- read.csv("data/TCRseq_metadata.csv") %>% select(UIN,sample)
div <- read.csv("data/Diversity_down-sampled_beta.csv") %>%
  filter(tissue == "TST_D7") %>%
  left_join(meta_tcr) %>%
  inner_join(meta) %>%
  mutate('Richness' = c(scale(richness,scale = T, center = T)),
         'Gini' = c(scale(gini, scale = T, center = T)),
         'Shannon' = c(scale(shannon, scale = T, center = T)),
         'Simpson' = c(scale(invsimpson, scale = T, center = T)),
         'Expanded' = c(scale(sum_expanded, scale = T, center = T))) %>%
  select(UIN,TB.exposure,Richness,Gini,Shannon,Simpson,Expanded) %>%
  pivot_longer(cols = c("Richness","Gini","Shannon","Simpson","Expanded"))
div$name <- factor(div$name, levels = c("Expanded","Gini","Richness","Shannon","Simpson"))

stats <- div %>%
  group_by(name) %>%
  wilcox_test(data = ., value~TB.exposure) %>%
  add_xy_position() # n=17 vs. 43

pC <- ggplot(div, aes(x=TB.exposure,y=value)) +
  geom_boxplot() +
  facet_grid(~name) +
  stat_pvalue_manual(stats) +
  My_Theme +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("") +
  ylab("Z score")
pC

## within sample co-incidence (TST_D7) ####
wsc <- read.csv("data/pc_withindonor_beta_down-sampled.csv") %>%
  select(sample,UIN,tissue,logpc) %>%
  filter(tissue == "TST_D7") %>%
  inner_join(meta)

stats <- wsc %>%
  wilcox_test(data = ., logpc~TB.exposure) %>%
  add_xy_position() # n=17 vs. 43

pD <- ggplot(wsc, aes(x=TB.exposure,y=logpc)) +
  geom_boxplot() +
  stat_pvalue_manual(stats) +
  My_Theme +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("") +
  ylab("Within donor coincidence")
pD

## cross sample co-incidence (TST_D7) ####
csc <- read.csv("data/pc_crossdonor_beta_down-sampled.csv") %>%
  filter(tissue == "TST_D7") %>%
  dplyr::rename(sample = "group1") %>%
  left_join(meta_tcr) %>%
  left_join(meta) %>%
  dplyr::rename(group1 = "sample", UIN1 = "UIN", TB.exposure1 = "TB.exposure") %>%
  dplyr::rename(sample = "group2") %>%
  left_join(meta_tcr) %>%
  left_join(meta) %>%
  dplyr::rename(group2 = "sample", UIN2 = "UIN", TB.exposure2 = "TB.exposure") %>%
  na.omit()

csc1 <- csc %>%
  select(logpc,TB.exposure1,TB.exposure2) %>%
  mutate(TB.exposure = ifelse(TB.exposure1 == TB.exposure2, TB.exposure1, NA)) %>%
  na.omit()

stats <- csc1 %>%  
  wilcox_test(data = ., logpc~TB.exposure) %>%
  add_xy_position() # n=136 vs. 903 pairwise comparisons

pE <- ggplot(csc1, aes(x=TB.exposure,y=logpc)) +
  geom_boxplot() +
  stat_pvalue_manual(stats) +
  My_Theme +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("") +
  ylab("Between donor coincidence")
pE

## abundance of Mtb-reactive TCRs in day 7 TSTs ####
### public Mtb-reactive CDR3s ####
a <- read.csv("data/Published-Ag-abundance_down-sampled_expanded_gr0_beta.csv") %>%
  filter(tissue == "TST_D7" & Antigen == "Mtb") %>%
  mutate(Clone.Size = "All")
b <- read.csv("data/Published-Ag-abundance_down-sampled_expanded_gr1_beta.csv") %>%
  filter(tissue == "TST_D7" & Antigen == "Mtb") %>%
  mutate(Clone.Size = "Expanded")
summary <- rbind(a,b) %>%
  left_join(meta_tcr) %>%
  inner_join(meta)

stats <- summary %>%
  group_by(Clone.Size) %>%
  wilcox_test(data = ., pct~TB.exposure) %>%
  add_xy_position() # n=17 vs. 43

pF <- ggplot(summary, aes(x=TB.exposure,y=pct)) +
  geom_boxplot() +
  facet_grid(~Clone.Size) +
  stat_pvalue_manual(stats) +
  My_Theme +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("") +
  ylab("% public Mtb-reactive CDR3s")
pF

### metaclones ####
a <- read.csv("data/summary_metaclone-abundance_down-sampled_beta_expanded_gr0.csv") %>%
  filter(tissue == "TST_D7") %>%
  mutate(Clone.Size = "All")
b <- read.csv("data/summary_metaclone-abundance_down-sampled_beta_expanded_gr1.csv") %>%
  filter(tissue == "TST_D7") %>%
  mutate(Clone.Size = "Expanded")

summary <- rbind(a,b) %>%
  left_join(meta_tcr) %>%
  inner_join(meta)

stats <- summary %>%
  group_by(Clone.Size) %>%
  wilcox_test(data = ., mc.pct~TB.exposure) %>%
  add_xy_position() # n=17 vs. 43

pG <- ggplot(summary, aes(x=TB.exposure,y=mc.pct)) +
  geom_boxplot() +
  facet_grid(~Clone.Size) +
  stat_pvalue_manual(stats) +
  My_Theme +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("") +
  ylab("% metaclone TCRs")
pG

### expanded TCRs ####
exp <- read.csv("data/summary_pct_expanded_tcrs_down-sampled_all-tissue.csv", row.names = 1) %>%
  filter(tissue == "TST_D7") %>%
  mutate(UIN = str_split_i(sample,"_",1)) %>%
  select(UIN,starts_with("pct")) %>%
  inner_join(meta) %>%
  pivot_longer(cols = starts_with("pct")) %>%
  mutate(name = recode(name,
                       pct_gr_1 = ">1",
                       pct_gr_2 = ">2",
                       pct_gr_3 = ">3",
                       pct_gr_4 = ">4"))

stats <- exp %>%
  group_by(name) %>%
  wilcox_test(data = ., value~TB.exposure) %>%
  add_xy_position() # n=17 vs. 43

pH <- ggplot(exp, aes(x=TB.exposure,y=value)) +
  geom_boxplot() +
  facet_grid(~name) +
  stat_pvalue_manual(stats) +
  My_Theme +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("") +
  ylab("% expanded TCRs")
pH

#### summary figure ####
r1 <- ggarrange(pA,pB,
                ncol = 2,
                widths = c(1.25,5.75),
                labels = c("A","B"),
                font.label = list(size=10, face = "bold", colour = "black"))
r2 <- ggarrange(pC,pD,pE,
                ncol = 3,
                widths = c(4.5,1.25,1.25),
                labels = c("C","D","E"),
                font.label = list(size = 10, face = "bold", colour = "black"))
r3 <- ggarrange(pF,pG,pH,
                ncol = 3,
                widths = c(2.25,2.25,3),
                labels = c("F","G","H"),
                font.label = list(size = 10, face = "bold", colour = "black"))

ggarrange(r1,r2,r3,
          nrow = 3,
          heights = c(1.1,0.95,0.95))

ggsave("figures/FigureS15.svg", 
       units = "cm", width = 17, height =22 , dpi=300)
