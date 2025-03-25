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
d <- read.csv("data/Diversity_down-sampled_beta.csv")

dat <- d %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST")) %>%
  select(richness,gini,shannon,invsimpson,sum_expanded,tissue)

# scaled data
zdat <- dat %>%
  mutate(richness = c(scale(richness,scale = T, center = T)),
         gini = c(scale(gini, scale = T, center = T)),
         shannon = c(scale(shannon, scale = T, center = T)),
         invsimpson = c(scale(invsimpson, scale = T, center = T)),
         sum_expanded = c(scale(sum_expanded, scale = T, center = T))) %>%
  ungroup()

# prep for plotting
zdat.long <- zdat %>%
  pivot_longer(cols = richness:sum_expanded,
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         richness = "Richness",
                         gini = "Gini index",
                         shannon = "Shannon diversity",
                         invsimpson = "Simpson diversity",
                         sum_expanded = "No. expanded TCRs"))

zdat.long$Metric <- factor(zdat.long$Metric, levels = c(
  "No. expanded TCRs",
  "Gini index",
  "Richness",
  "Shannon diversity",
  "Simpson diversity"))

# stats
stats <- zdat.long %>%
  group_by(Metric) %>%
  pairwise_wilcox_test(data = ., Value~tissue, p.adjust.method = "fdr")

# plot (save as svg 1000x350)
ggplot(zdat.long,aes(tissue,Value))+
  geom_boxplot(colour="blue",outliers = F)+
  facet_wrap(~Metric,ncol=5)+
  labs(x = "Sample", 
       y = "Z score")+
  My_Theme + 
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  stat_pvalue_manual(stats, label = "p.adj.signif", y.position = c(3.5,4,4.5))
