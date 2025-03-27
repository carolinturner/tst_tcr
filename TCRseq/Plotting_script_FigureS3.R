library(tidyverse)
library(ggpubr)
library(rstatix)

#My_Theme
t = 8 #size of text
m = 4 #size of margin around text
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
  legend.position = "right", legend.justification = "top",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0,0,0,0),
  legend.box.spacing = unit(c(0,0,0,0.1),"cm"),
  plot.margin = unit(c(0.5,0.6,0,0.2),"cm")
)

# Figure S3A: down-sampled, alpha ####
d <- read.csv("data/Diversity_down-sampled_alpha.csv") %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST")) %>%
  select(richness,gini,shannon,invsimpson,sum_expanded,tissue)

# scaled data
zdat <- d %>%
  mutate('Richness' = c(scale(richness,scale = T, center = T)),
         'Gini index' = c(scale(gini, scale = T, center = T)),
         'Shannon diversity' = c(scale(shannon, scale = T, center = T)),
         'Simpson diversity' = c(scale(invsimpson, scale = T, center = T)),
         'No. expanded TCRs' = c(scale(sum_expanded, scale = T, center = T)))

# prep for plotting
zdat.long <- zdat %>%
  pivot_longer(cols = 'Richness':'No. expanded TCRs',
               names_to = "Metric",
               values_to = "Value") 

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

# plot 
pS3A <- ggplot(zdat.long,aes(tissue,Value))+
  geom_boxplot(colour="red",outliers = F)+
  facet_wrap(~Metric,ncol=5)+
  labs(x = "Sample", 
       y = "Z score")+
  My_Theme + 
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  stat_pvalue_manual(stats, label = "p.adj.signif", y.position = c(3.5,4.5,5.5), size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  ggtitle("alpha (down-sampled)")

# Figure S3B: full repertoires, alpha ####
d <- read.csv("data/Diversity_full-repertoires_alpha.csv") %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST")) %>%
  select(richness,gini,shannon,invsimpson,sum_expanded,tissue)

# scaled data
zdat <- d %>%
  mutate('Richness' = c(scale(richness,scale = T, center = T)),
         'Gini index' = c(scale(gini, scale = T, center = T)),
         'Shannon diversity' = c(scale(shannon, scale = T, center = T)),
         'Simpson diversity' = c(scale(invsimpson, scale = T, center = T)),
         'No. expanded TCRs' = c(scale(sum_expanded, scale = T, center = T))) 

# prep for plotting
zdat.long <- zdat %>%
  pivot_longer(cols = 'Richness':'No. expanded TCRs',
               names_to = "Metric",
               values_to = "Value") 

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

# plot
pS3B <- ggplot(zdat.long,aes(tissue,Value))+
  geom_boxplot(colour="red",outliers = F)+
  facet_wrap(~Metric,ncol=5)+
  labs(x = "Sample", 
       y = "Z score")+
  My_Theme + 
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  stat_pvalue_manual(stats, label = "p.adj.signif", y.position = c(5,6,7), size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  ggtitle("alpha (full repertoires)")

# Figure S3C: full repertoires, beta ####
d <- read.csv("data/Diversity_full-repertoires_beta.csv") %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST")) %>%
  select(richness,gini,shannon,invsimpson,sum_expanded,tissue)

# scaled data
zdat <- d %>%
  mutate('Richness' = c(scale(richness,scale = T, center = T)),
         'Gini index' = c(scale(gini, scale = T, center = T)),
         'Shannon diversity' = c(scale(shannon, scale = T, center = T)),
         'Simpson diversity' = c(scale(invsimpson, scale = T, center = T)),
         'No. expanded TCRs' = c(scale(sum_expanded, scale = T, center = T)))

# prep for plotting
zdat.long <- zdat %>%
  pivot_longer(cols = 'Richness':'No. expanded TCRs',
               names_to = "Metric",
               values_to = "Value")

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

# plot
pS3C <- ggplot(zdat.long,aes(tissue,Value))+
  geom_boxplot(colour="blue",outliers = F)+
  facet_wrap(~Metric,ncol=5)+
  labs(x = "Sample", 
       y = "Z score")+
  My_Theme + 
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  stat_pvalue_manual(stats, label = "p.adj.signif", y.position = c(5,6,7), size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  ggtitle("beta (full repertoires)")

# assemble figure ####
ggarrange(pS3A,pS3B,pS3C,
          nrow=3,
          labels = list("A","B","C"),
          font.label = list(size = 10, face = "bold", colour = "black"))
ggsave("FigureS3.svg", 
       units = "cm", width = 17, height =17, dpi=300)
