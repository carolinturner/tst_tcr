library(tidyverse)
library(rstatix)
library(ggpubr)

#My_Theme
t = 8 #size of text
m = 2 #size of margin around text
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
  plot.margin = unit(c(0.5,0.2,0,0.1),"cm"),
  panel.spacing = unit(0.05,"cm")
)


# Figure S4A: down-sampled, alpha ####
a <- read.csv("data/Published-Ag-abundance_down-sampled_expanded_gr0_alpha.csv") %>% mutate(Clone.Size = "All TCRs")
b <- read.csv("data/Published-Ag-abundance_down-sampled_expanded_gr1_alpha.csv") %>% mutate(Clone.Size = "Expanded TCRs")

# number of published sequences
ref <- read.csv("data/FileS3.csv") %>% filter(chain == "alpha")
CMV <- length(ref %>% filter(reactivity == "CMV") %>% pull(CDR3) %>% unique())
CMV <- format(c(CMV),big.mark=",", trim=TRUE)
EBV <- length(ref %>% filter(reactivity == "EBV") %>% pull(CDR3) %>% unique())
EBV <- format(c(EBV),big.mark=",", trim=TRUE)
Mtb <- length(ref %>% filter(reactivity == "Mtb") %>% pull(CDR3) %>% unique())
Mtb <- format(c(Mtb),big.mark=",", trim=TRUE)

summary <- rbind(a,b) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         Antigen = recode(Antigen,
                          CMV = paste0("CMV\n(",CMV," CDR3s)"),
                          EBV = paste0("EBV\n(",EBV," CDR3s)"),
                          Mtb = paste0("Mtb\n(",Mtb," CDR3s)")))
summary$Antigen <- factor(summary$Antigen,
                          levels = c(paste0("Mtb\n(",Mtb," CDR3s)"),
                                     paste0("CMV\n(",CMV," CDR3s)"),
                                     paste0("EBV\n(",EBV," CDR3s)")))

# stats
stats <- summary %>%
  group_by(Antigen,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot 
pS4A <- ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="red")+
  facet_grid(Clone.Size~Antigen)+
  scale_y_log10()+
  labs(x="Sample",
       y="% of all CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank())+
  expand_limits(y=130)+
  stat_pvalue_manual(stats,
                     label = "p.adj.signif",
                     y.position = c(log10(45),log10(70),log10(110)),
                     size = 2)

# Figure S4B: full repertoires, alpha ####
a <- read.csv("data/Published-Ag-abundance_full-repertoires_expanded_gr0_alpha.csv") %>% mutate(Clone.Size = "All TCRs")
b <- read.csv("data/Published-Ag-abundance_full-repertoires_expanded_gr1_alpha.csv") %>% mutate(Clone.Size = "Expanded TCRs")

# number of published sequences
ref <- read.csv("data/FileS3.csv") %>% filter(chain == "alpha")
CMV <- length(ref %>% filter(reactivity == "CMV") %>% pull(CDR3) %>% unique())
CMV <- format(c(CMV),big.mark=",", trim=TRUE)
EBV <- length(ref %>% filter(reactivity == "EBV") %>% pull(CDR3) %>% unique())
EBV <- format(c(EBV),big.mark=",", trim=TRUE)
Mtb <- length(ref %>% filter(reactivity == "Mtb") %>% pull(CDR3) %>% unique())
Mtb <- format(c(Mtb),big.mark=",", trim=TRUE)

summary <- rbind(a,b) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         Antigen = recode(Antigen,
                          CMV = paste0("CMV\n(",CMV," CDR3s)"),
                          EBV = paste0("EBV\n(",EBV," CDR3s)"),
                          Mtb = paste0("Mtb\n(",Mtb," CDR3s)")))
summary$Antigen <- factor(summary$Antigen,
                          levels = c(paste0("Mtb\n(",Mtb," CDR3s)"),
                                     paste0("CMV\n(",CMV," CDR3s)"),
                                     paste0("EBV\n(",EBV," CDR3s)")))

# stats
stats <- summary %>%
  group_by(Antigen,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot 
pS4B <- ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="red")+
  facet_grid(Clone.Size~Antigen)+
  scale_y_log10()+
  labs(x="Sample",
       y="% of all CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank())+
  expand_limits(y=130)+
  stat_pvalue_manual(stats,
                     label = "p.adj.signif",
                     y.position = c(log10(45),log10(70),log10(110)),
                     size=2)

# Figure S4C: full repertoires, beta ####
a <- read.csv("data/Published-Ag-abundance_full-repertoires_expanded_gr0_beta.csv") %>% mutate(Clone.Size = "All TCRs")
b <- read.csv("data/Published-Ag-abundance_full-repertoires_expanded_gr1_beta.csv") %>% mutate(Clone.Size = "Expanded TCRs")

# number of published sequences
ref <- read.csv("data/FileS3.csv") %>% filter(chain == "beta")
CMV <- length(ref %>% filter(reactivity == "CMV") %>% pull(CDR3) %>% unique())
CMV <- format(c(CMV),big.mark=",", trim=TRUE)
EBV <- length(ref %>% filter(reactivity == "EBV") %>% pull(CDR3) %>% unique())
EBV <- format(c(EBV),big.mark=",", trim=TRUE)
Mtb <- length(ref %>% filter(reactivity == "Mtb") %>% pull(CDR3) %>% unique())
Mtb <- format(c(Mtb),big.mark=",", trim=TRUE)

summary <- rbind(a,b) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"),
         Antigen = recode(Antigen,
                          CMV = paste0("CMV\n(",CMV," CDR3s)"),
                          EBV = paste0("EBV\n(",EBV," CDR3s)"),
                          Mtb = paste0("Mtb\n(",Mtb," CDR3s)")))
summary$Antigen <- factor(summary$Antigen,
                          levels = c(paste0("Mtb\n(",Mtb," CDR3s)"),
                                     paste0("CMV\n(",CMV," CDR3s)"),
                                     paste0("EBV\n(",EBV," CDR3s)")))

# stats
stats <- summary %>%
  group_by(Antigen,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot 
pS4C <- ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="blue")+
  facet_grid(Clone.Size~Antigen)+
  scale_y_log10(limits = c(0.016,100))+
  labs(x="Sample",
       y="% of all CDR3s") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats,
                     label = "p.adj.signif",
                     y.position = c(log10(20),log10(40),log10(80)),
                     size = 2)

# assemble figure ####
r1 <- ggarrange(pS4A,
                nrow = 1,
                ncol = 2,
                labels = list("A"),
                font.label = list(size = 10, face = "bold", colour = "black"))
r2 <- ggarrange(pS4B,pS4C,
                nrow=1,
                labels = list("B","C"),
                font.label = list(size = 10, face = "bold", colour = "black"))
ggarrange(r1,r2,
          nrow=2,
          font.label = list(size = 10, face = "bold", colour = "black"))
ggsave("figures/FigureS4.svg", 
       units = "cm", width = 17, height = 18 , dpi=300)
