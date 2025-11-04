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
  panel.spacing = unit(0.05,"cm"),
  plot.margin = unit(c(0.5,0.6,0,0.2),"cm")
)

summary <- read.csv("data/invariant_alpha_tcrs.csv")
df <- summary %>%
  filter(tissue %in% c("Blood","TST_D2","TST_D7")) %>%
  mutate(dataset = recode(dataset,
                          "down" = "down-sampled",
                          "full" = "full repertoires"))
df$type <- factor(df$type, levels = c("MAIT","iNKT","GEM"))

# stats
stats <- df %>%
  group_by(dataset,type) %>%
  pairwise_wilcox_test(data = ., inv.pct~tissue, p.adjust.method = "fdr") %>%
  add_xy_position()

# plot 
p <- ggplot(df, aes(x=tissue,y=inv.pct))+
  geom_boxplot(colour="red")+
  facet_grid(type~dataset, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x="Sample",
       y="% of all TCRs") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats, label = "p.adj.signif",size = 3)
p

ggsave("figures/FigureS5.svg", 
       units = "cm", width = 10, height =10 , dpi=300)
