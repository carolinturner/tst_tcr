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

# Figure S3A: down-sampled, alpha ####
a <- read.csv("data/Published-Ag-abundance_down-sampled_expanded_gr0_alpha.csv")
b <- read.csv("data/Published-Ag-abundance_down-sampled_expanded_gr1_alpha.csv")

a <- a %>% mutate(Clone.Size = "All CDR3s")
b <- b %>% mutate(Clone.Size = "Expanded CDR3s")

# number of published sequences
ref <- read.csv("data/TableS2.csv") %>%
  filter(chain == "alpha")
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

# stats
stats <- summary %>%
  group_by(Antigen,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot (save as svg 600x500)
ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="red")+
  facet_grid(Clone.Size~Antigen)+
  scale_y_log10()+
  labs(x="Sample",
       y="% of all CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(20),log10(40),log10(80)))

# Figure S3B: full repertoires, beta ####
a <- read.csv("data/Published-Ag-abundance_full-repertoires_expanded_gr0_beta.csv")
b <- read.csv("data/Published-Ag-abundance_full-repertoires_expanded_gr1_beta.csv")

a <- a %>% mutate(Clone.Size = "All CDR3s")
b <- b %>% mutate(Clone.Size = "Expanded CDR3s")

# number of published sequences
ref <- read.csv("data/TableS2.csv") %>%
  filter(chain == "beta")
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

# stats
stats <- summary %>%
  group_by(Antigen,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot (save as svg 600x500)
ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="blue")+
  facet_grid(Clone.Size~Antigen)+
  scale_y_log10(limits = c(0.016,100))+
  labs(x="Sample",
       y="% of all CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(20),log10(40),log10(80)))

# Figure S3C: full repertoires, alpha ####
a <- read.csv("data/Published-Ag-abundance_full-repertoires_expanded_gr0_alpha.csv")
b <- read.csv("data/Published-Ag-abundance_full-repertoires_expanded_gr1_alpha.csv")

a <- a %>% mutate(Clone.Size = "All CDR3s")
b <- b %>% mutate(Clone.Size = "Expanded CDR3s")

# number of published sequences
ref <- read.csv("data/TableS2.csv") %>%
  filter(chain == "alpha")
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

# stats
stats <- summary %>%
  group_by(Antigen,Clone.Size) %>%
  pairwise_wilcox_test(data = ., pct~tissue, p.adjust.method = "fdr") 

# plot (save as svg 600x500)
ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="red")+
  facet_grid(Clone.Size~Antigen)+
  scale_y_log10()+
  labs(x="Sample",
       y="% of all CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(20),log10(40),log10(80)))
