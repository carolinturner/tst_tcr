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
  panel.background = element_rect(fill = "gray97")
)

# Figure S9A: down-sampled beta ###
a <- read.csv("data/Summary_percentage_beta_size-matched-mc-mtb_search_down-sampled_expanded_gr0.csv")
b <- read.csv("data/Summary_percentage_beta_size-matched-mc-mtb_search_down-sampled_expanded_gr1.csv")

a <- a %>% mutate(Clone.Size = "All CDR3s")
b <- b %>% mutate(Clone.Size = "Expanded CDR3s")

summary <- rbind(a,b) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         TCR = recode(TCR,
                      published = "published (1,392 CDR3s)",
                      metaclone = "metaclone (1,392 CDR3s)"))

# stats
stats <- summary %>%
  group_by(TCR,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot (save as svg 500x500)
ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="blue")+
  facet_grid(Clone.Size~TCR)+
  scale_y_log10(expand = expansion(mult=c(0,0.1)))+
  labs(x="Sample",
       y="% of all CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(10),log10(20),log10(40)))


# Figure S9B: down-sampled: alpha ####
a <- read.csv("data/Summary_percentage_alpha_size-matched-mc-mtb_search_down-sampled_expanded_gr0.csv")
b <- read.csv("data/Summary_percentage_alpha_size-matched-mc-mtb_search_down-sampled_expanded_gr1.csv")

a <- a %>% mutate(Clone.Size = "All CDR3s")
b <- b %>% mutate(Clone.Size = "Expanded CDR3s")

summary <- rbind(a,b) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         TCR = recode(TCR,
                      published = "published (545 CDR3s)",
                      metaclone = "metaclone (545 CDR3s)"))

# stats
stats <- summary %>%
  group_by(TCR,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot (save 500x500)
ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="red")+
  facet_grid(Clone.Size~TCR)+
  scale_y_log10(expand = expansion(mult=c(0,0.1)))+
  labs(x="Sample",
       y="% of all CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(10),log10(20),log10(40)))

# Figure S9C: full repertoires beta ####
a <- read.csv("data/Summary_percentage_beta_size-matched-mc-mtb_search_full-repertoires_expanded_gr0.csv")
b <- read.csv("data/Summary_percentage_beta_size-matched-mc-mtb_search_full-repertoires_expanded_gr1.csv")

a <- a %>% mutate(Clone.Size = "All CDR3s")
b <- b %>% mutate(Clone.Size = "Expanded CDR3s")

summary <- rbind(a,b) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         TCR = recode(TCR,
                      published = "published (1,392 CDR3s)",
                      metaclone = "metaclone (1,392 CDR3s)"))

# stats
stats <- summary %>%
  group_by(TCR,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot (save as svg 500x500)
ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="blue")+
  facet_grid(Clone.Size~TCR)+
  scale_y_log10(expand = expansion(mult=c(0,0.1)))+
  labs(x="Sample",
       y="% of all CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(10),log10(20),log10(40)))

# Figure S9D: full repertoires alpha ###
a <- read.csv("data/Summary_percentage_alpha_size-matched-mc-mtb_search_full-repertoires_expanded_gr0.csv")
b <- read.csv("data/Summary_percentage_alpha_size-matched-mc-mtb_search_full-repertoires_expanded_gr1.csv")

a <- a %>% mutate(Clone.Size = "All CDR3s")
b <- b %>% mutate(Clone.Size = "Expanded CDR3s")

summary <- rbind(a,b) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         TCR = recode(TCR,
                      published = "published (545 CDR3s)",
                      metaclone = "metaclone (545 CDR3s)"))

# stats
stats <- summary %>%
  group_by(TCR,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot (save as svg 500x500)
ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="red")+
  facet_grid(Clone.Size~TCR)+
  scale_y_log10(expand = expansion(mult=c(0,0.1)))+
  labs(x="Sample",
       y="% of all CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(10),log10(20),log10(40)))
