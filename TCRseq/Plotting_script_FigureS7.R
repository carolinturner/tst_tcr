library(tidyverse)
library(ggplotify)
library(pheatmap)
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
  plot.margin = unit(c(0.5,0.2,0,0.2),"cm"),
  panel.spacing.x = unit(0.05,"cm")
)

# Figure S7A ####
dat.a <- read.csv("data/Publicness_in-vitro_PPD.csv",row.names = 1) %>%
  select(-beta)

# plot 
pS7A <- ggplot(dat.a, aes(x=publicness,y=alpha))+
  geom_bar(stat = "identity",fill="red",colour="red",width=0.8)+
  scale_x_continuous(n.breaks = 11, limits = c(0.6,12))+
  scale_y_log10()+
  labs(x = "Number of participants", 
       y = "No. unique CDR3s") +
  My_Theme

# Figure S7B ####
dat <- read.csv("data/Heatmap_data_alpha.csv",row.names = 1)

hm <- as.ggplot(
  pheatmap(dat,
           fontsize = 10,
           color=colorRampPalette(c("grey90", "red"))(2),
           legend = F,
           border_color = "black",
           cluster_rows = T,
           cluster_cols = T,
           clustering_method = "ward.D2",
           show_rownames = F,
           show_colnames = F,
           treeheight_row = 20,
           treeheight_col = 20))
pS7B <- ggarrange(hm) %>%
  annotate_figure(bottom = text_grob("Participants", size = 8, face = "bold"),
                  left = text_grob("CDR3 sequences", size = 8, face = "bold", rot = 90))


# Figure S7C-D: down-sampled repertoires ####
a <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr0.csv")
b <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr1.csv")
c <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr2.csv")
d <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr3.csv")
e <- read.csv("data/Summary_down-sampled_private-Ag-abundance_expanded_gr4.csv")

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "TST-2",
                         TST_D7 = "TST-7"),
         Antigen = recode(Antigen,
                          PPD = "PPD-react.",
                          TT = "TT-react."))

# alpha chain data
alpha <- summary %>% filter(chain == "alpha")

blood.stats <- alpha %>%
  filter(tissue == "Blood" & Antigen == "PPD-react.")
blood.pct.all <- blood.stats %>%
  pairwise_wilcox_test(data = ., pct~Clone.Size, paired = TRUE, p.adjust.method = "fdr")
blood.pct.uniq <- blood.stats %>%
  pairwise_wilcox_test(data = ., pct.unique~Clone.Size, paired = TRUE, p.adjust.method = "fdr")

tstd2.stats <- alpha %>%
  filter(tissue == "TST-2" & Antigen == "PPD-react.")
tstd2.pct.all <- tstd2.stats %>%
  pairwise_wilcox_test(data = ., pct~Clone.Size, paired = TRUE, p.adjust.method = "fdr")
tstd2.pct.uniq <- tstd2.stats %>%
  pairwise_wilcox_test(data = ., pct.unique~Clone.Size, paired = TRUE, p.adjust.method = "fdr")

# Figure S7C: all CDR3s
pS7C <- ggplot(alpha, aes(x=tissue,y=pct))+
  geom_boxplot(colour="red")+
  facet_grid(Antigen~Clone.Size)+
  ylim(0,100)+
  labs(y="% all CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  ggtitle("alpha (down-sampled)")

# Figure S7D: unique CDR3s
pS7D <- ggplot(alpha, aes(x=tissue,y=pct.unique))+
  geom_boxplot(colour="red")+
  facet_grid(Antigen~Clone.Size)+
  ylim(0,100)+
  labs(y="% unique CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  ggtitle("alpha (down-sampled)")

# Figure S7E-H: full repertoires ####
a <- read.csv("data/Summary_full-repertoires_private-Ag-abundance_expanded_gr0.csv")
b <- read.csv("data/Summary_full-repertoires_private-Ag-abundance_expanded_gr1.csv")
c <- read.csv("data/Summary_full-repertoires_private-Ag-abundance_expanded_gr2.csv")
d <- read.csv("data/Summary_full-repertoires_private-Ag-abundance_expanded_gr3.csv")
e <- read.csv("data/Summary_full-repertoires_private-Ag-abundance_expanded_gr4.csv")

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "TST-2",
                         TST_D7 = "TST-7"),,
         Antigen = recode(Antigen,
                          PPD = "PPD-react.",
                          TT = "TT-react."))

# split by chain
beta <- summary %>% filter(chain == "beta")
blood.stats <- beta %>%
  filter(tissue == "Blood" & Antigen == "PPD-react.")
blood.pct.all <- blood.stats %>%
  pairwise_wilcox_test(data = ., pct~Clone.Size, paired = TRUE, p.adjust.method = "fdr")
blood.pct.uniq <- blood.stats %>%
  pairwise_wilcox_test(data = ., pct.unique~Clone.Size, paired = TRUE, p.adjust.method = "fdr")

tstd2.stats <- beta %>%
  filter(tissue == "TST-2" & Antigen == "PPD-react.")
tstd2.pct.all <- tstd2.stats %>%
  pairwise_wilcox_test(data = ., pct~Clone.Size, paired = TRUE, p.adjust.method = "fdr")
tstd2.pct.uniq <- tstd2.stats %>%
  pairwise_wilcox_test(data = ., pct.unique~Clone.Size, paired = TRUE, p.adjust.method = "fdr")

alpha <- summary %>% filter(chain == "alpha")
blood.stats <- alpha %>%
  filter(tissue == "Blood" & Antigen == "PPD-react.")
blood.pct.all <- blood.stats %>%
  pairwise_wilcox_test(data = ., pct~Clone.Size, paired = TRUE, p.adjust.method = "fdr")
blood.pct.uniq <- blood.stats %>%
  pairwise_wilcox_test(data = ., pct.unique~Clone.Size, paired = TRUE, p.adjust.method = "fdr")

tstd2.stats <- alpha %>%
  filter(tissue == "TST-2" & Antigen == "PPD-react.")
tstd2.pct.all <- tstd2.stats %>%
  pairwise_wilcox_test(data = ., pct~Clone.Size, paired = TRUE, p.adjust.method = "fdr")
tstd2.pct.uniq <- tstd2.stats %>%
  pairwise_wilcox_test(data = ., pct.unique~Clone.Size, paired = TRUE, p.adjust.method = "fdr")


# Figure S7E: all CDR3s alpha
pS7E <- ggplot(alpha, aes(x=tissue,y=pct))+
  geom_boxplot(colour="red")+
  facet_grid(Antigen~Clone.Size)+
  ylim(0,100)+
  labs(y="%all CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  ggtitle("alpha (full repertoire)")

# Figure S7F: unique CDR3s alpha
pS7F <- ggplot(alpha, aes(x=tissue,y=pct.unique))+
  geom_boxplot(colour="red")+
  facet_grid(Antigen~Clone.Size)+
  ylim(0,100)+
  labs(y="% unique CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  ggtitle("alpha (full repertoire)")

# Figure S7G: all CDR3s beta
pS7G <- ggplot(beta, aes(x=tissue,y=pct))+
  geom_boxplot(colour="blue")+
  facet_grid(Antigen~Clone.Size)+
  ylim(0,100)+
  labs(y="% all CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  ggtitle("beta (full repertoire)")

# Figure S7H: unique CDR3s beta
pS7H <- ggplot(beta, aes(x=tissue,y=pct.unique))+
  geom_boxplot(colour="blue")+
  facet_grid(Antigen~Clone.Size)+
  ylim(0,100)+
  labs(y="% unique CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  ggtitle("beta (full repertoire)")

# assemble figure ####
r1 <- ggarrange(pS7A,pS7B,
                ncol=3,
                widths = c(1.3,1,0.7),
                labels = list("A","B"),
                font.label = list(size = 10, face = "bold", colour = "black"))
r2 <- ggarrange(pS7C,pS7D,pS7E,pS7F,pS7G,pS7H,
                ncol=2,nrow=3,
                labels = list("C","D","E","F","G","H"),
                font.label = list(size = 10, face = "bold", colour = "black"))
ggarrange(r1,r2,
          nrow = 2,
          heights = c(0.7,3.5))
ggsave("figures/FigureS7.png", 
       units = "cm", width = 17, height =25 , dpi=300)
