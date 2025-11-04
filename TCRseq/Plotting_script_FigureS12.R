library(tidyverse)
library(rstatix)
library(ggpubr)

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
  panel.spacing = unit(0.05,"cm"),
  plot.margin = unit(c(0.2,0.6,0.5,0.5),"cm")
)

# Figure S12A: down-sampled beta ###
a <- read.csv("data/Summary_percentage_beta_size-matched-mc-mtb_search_down-sampled_expanded_gr0.csv")
b <- read.csv("data/Summary_percentage_beta_size-matched-mc-mtb_search_down-sampled_expanded_gr1.csv")

a <- a %>% mutate(Clone.Size = "All TCRs")
b <- b %>% mutate(Clone.Size = "Expanded TCRs")

summary <- rbind(a,b) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         TCR = recode(TCR,
                      published = "published (1392)",
                      metaclone = "metaclone (1392)"))

# stats
stats <- summary %>%
  group_by(TCR,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot 
pS12A <- ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="blue")+
  facet_grid(Clone.Size~TCR)+
  scale_y_log10(expand = expansion(mult=c(0,0.1)))+
  labs(x="Sample",
       y="% all TCRs") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(10),log10(20),log10(40)), size = 2)+
  ggtitle("beta (down-sampled)")


# Figure S12B: full repertoires beta ####
a <- read.csv("data/Summary_percentage_beta_size-matched-mc-mtb_search_full-repertoires_expanded_gr0.csv")
b <- read.csv("data/Summary_percentage_beta_size-matched-mc-mtb_search_full-repertoires_expanded_gr1.csv")

a <- a %>% mutate(Clone.Size = "All TCRs")
b <- b %>% mutate(Clone.Size = "Expanded TCRs")

summary <- rbind(a,b) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         TCR = recode(TCR,
                      published = "published (1392)",
                      metaclone = "metaclone (1392)"))

# stats
stats <- summary %>%
  group_by(TCR,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot 
pS12B <- ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="blue")+
  facet_grid(Clone.Size~TCR)+
  scale_y_log10(expand = expansion(mult=c(0,0.1)))+
  labs(x="Sample",
       y="% all TCRs") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(10),log10(20),log10(40)), size = 2)+
  ggtitle("beta (full repertoire)")

# assemble figure ####
ggarrange(pS12A,pS12B,
          ncol = 2,
          labels = c("A","B"),
          font.label = list(size = 10, face = "bold", colour = "black"))
ggsave("figures/FigureS12.svg", 
       units = "cm", width = 17, height =9 , dpi=300)
