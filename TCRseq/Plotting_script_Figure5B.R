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

a <- read.csv("data/summary_metaclone-abundance_down-sampled_expanded_gr0.csv")
b <- read.csv("data/summary_metaclone-abundance_down-sampled_expanded_gr1.csv")

# beta chain data
c <- a %>% filter(Chain == "beta") %>% na.omit() %>% mutate(Clone.Size = "All CDR3s")
d <- b %>% filter(Chain == "beta") %>% na.omit() %>% mutate(Clone.Size = "Expanded CDR3s")

summary <- rbind(c,d) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"))
# stats
stats.all <- summary %>%
  group_by(Clone.Size) %>%
  pairwise_wilcox_test(data = .,mc.pct~tissue, p.adjust.method = "fdr")

# plot (save as svg 350x300)
ggplot(summary, aes(x=tissue,y=mc.pct))+
  geom_boxplot(colour="blue")+
  facet_wrap(~Clone.Size,ncol = 5)+
  scale_y_log10(limits = c(0.03,100))+
  labs(x="Sample",
       y="% of all CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  stat_pvalue_manual(stats.all, label = "p.adj.signif",y.position=c(log10(15),log10(30),log10(60)))
