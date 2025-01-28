library(tidyverse)
library(pheatmap)

#My_Theme
t = 10 #size of text
m = 10 #size of margin around text
tc = "black" #colour of text
My_Theme = theme(
  axis.title.x = element_text(size = t, face = "bold", margin = margin(t = m)),
  axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = 0, hjust = 0),
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
  legend.position = "right", legend.justification = "top"
)

# Figure 3A ####
dat <- read.csv("data/Publicness_in-vitro_PPD.csv",row.names = 1)
dat.b <- dat %>% select(-alpha)

# plot (save as svg 350x300)
ggplot(dat.b, aes(x=publicness,y=beta))+
  geom_bar(stat = "identity",fill="blue",colour="blue",width=0.8)+
  scale_x_continuous(n.breaks = 11, limits = c(0.6,12))+
  scale_y_log10()+
  labs(x = "Number of participants", 
       y = "Number of unique CDR3s") +
  My_Theme

# Figure 3B ####
dat <- read.csv("data/Heatmap_data_beta.csv",row.names = 1)

# plot
png("Figure3B.png",width=1808,height=1617,res=600)
pheatmap(dat,
         fontsize = 10,
         color=colorRampPalette(c("grey90", "blue"))(2),
         legend = F,
         border_color = "black",
         cluster_rows = T,
         cluster_cols = T,
         clustering_method = "ward.D2",
         show_rownames = F,
         show_colnames = F,
         treeheight_row = 20,
         treeheight_col = 20)
dev.off()

# Figure 3C-D ####
a <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr0.csv")
b <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr1.csv")
c <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr2.csv")
d <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr3.csv")
e <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr4.csv")


a <- a %>% mutate(Clone.Size = "Clone size >0")
b <- b %>% mutate(Clone.Size = "Clone size >1")
c <- c %>% mutate(Clone.Size = "Clone size >2")
d <- d %>% mutate(Clone.Size = "Clone size >3")
e <- e %>% mutate(Clone.Size = "Clone size >4")


summary <- rbind(a,b,c,d,e) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         Antigen = recode(Antigen,
                          PPD = "PPD-reactive",
                          TT = "TT-reactive"))

# beta chain data
beta <- summary %>% filter(chain == "beta")

# Figure 3C: all CDR3s (save as svg 800x400)
ggplot(beta, aes(x=tissue,y=pct))+
  geom_boxplot(colour="blue")+
  facet_grid(Antigen~Clone.Size)+
  ylim(0,100)+
  labs(x="Sample",
       y="% of all CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))


# Figure 3D: unique CDR3s (save as svg 800x400)
ggplot(beta, aes(x=tissue,y=pct.unique))+
  geom_boxplot(colour="blue")+
  facet_grid(Antigen~Clone.Size)+
  ylim(0,100)+
  labs(x="Sample",
       y="% of unique CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))

