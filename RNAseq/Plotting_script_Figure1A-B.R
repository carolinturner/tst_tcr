library(tidyverse)
library(ggpubr)
library(rstatix)

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

# data
d <- read.csv("data/Modules_Z-scores.csv")

d <- d %>%
  mutate(Stimulant=recode(Stimulant,
                          saline= "Saline",
                          TST_D2= "Day 2 TST",
                          TST_D7= "Day 7 TST"))

d$Stimulant <- factor(d$Stimulant, levels = c(
  "Saline",
  "Day 2 TST",
  "Day 7 TST"))
  
### plot proliferation module
mod <- "CCND1"

dat.ss <- d %>% 
  filter(module_name %in% mod) 

stats <- dat.ss %>%
  wilcox_test(data = ., module_Zscore~Stimulant, ref.group="Day 2 TST", p.adjust.method = "fdr") %>%
  add_xy_position()
  
ggplot(dat.ss,aes(Stimulant, module_Zscore))+
  geom_jitter(size = 1, alpha = 1, width = 0.1, height = 0)+
  geom_violin(trim = F, fill= "grey", alpha =0.5, colour = "red",  linewidth = 0.5)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Stimulation", 
       y = "Proliferation\nModule Z score") +
  My_Theme +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  stat_pvalue_manual(stats, label = "p.adj.signif")


### plot cell type modules
mod <- c("T_cell","CD4_T","CD8_T","NK","myeloid")

dat.ss <- d %>% 
  filter(module_name %in% mod) 

dat.ss$module_name <- factor(dat.ss$module_name, levels = c(
  "T_cell",
  "CD4_T",
  "CD8_T",
  "NK",
  "myeloid"))

stats <- dat.ss %>%
  group_by(module_name) %>%
  wilcox_test(data = ., module_Zscore~Stimulant,ref.group="Day 2 TST",p.adjust.method = "fdr") %>%
  add_xy_position()

ggplot(dat.ss,aes(Stimulant, module_Zscore))+
  geom_jitter(size = 0.5, alpha = 1, width = 0.1, height = 0)+
  geom_violin(trim = F, fill= "grey", alpha =0.5, colour = "red",  linewidth = 0.5)+
  facet_wrap(~module_name,
             labeller = as_labeller(c(
               T_cell="All T cells",
               CD4_T="CD4 T cells",
               CD8_T="CD8 T cells",
               NK="NK cells",
               myeloid="Myeloid cells")))+
  labs(x = "Stimulation", 
       y = "Module Z score") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  stat_pvalue_manual(stats, label = "p.adj.signif")
