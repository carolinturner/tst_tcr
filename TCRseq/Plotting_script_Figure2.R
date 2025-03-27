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
  panel.spacing = unit(0.05,"cm"),
  plot.margin = unit(c(0.5,0.6,0,0.2),"cm")
)



# Figure 2A: Diversity ####
# data
d <- read.csv("data/Diversity_down-sampled_beta.csv") %>%
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
p2A <- ggplot(zdat.long,aes(tissue,Value))+
  geom_boxplot(colour="blue",outliers = F)+
  facet_wrap(~Metric,ncol=5)+
  labs(x = "Sample", 
       y = "Z score")+
  My_Theme + 
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank()) +
  stat_pvalue_manual(stats, label = "p.adj.signif", y.position = c(3.8,4.5,5.2),size = 3)

# Figure 2B: published antigen reactivity ####
# data
a <- read.csv("data/Published-Ag-abundance_down-sampled_expanded_gr0_beta.csv") %>% mutate(Clone.Size = "All TCRs")
b <- read.csv("data/Published-Ag-abundance_down-sampled_expanded_gr1_beta.csv") %>% mutate(Clone.Size = "Expanded TCRs")

# number of published sequences
ref <- read.csv("data/TableS2.csv") %>% filter(chain == "beta")
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
p2B <- ggplot(summary, aes(x=tissue,y=pct))+
  geom_boxplot(colour="blue")+
  facet_grid(Clone.Size~Antigen)+
  scale_y_log10(limits = c(0.016,100))+
  labs(x="Sample",
       y="% of all CDR3s (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats, label = "p.adj.signif",y.position = c(log10(10),log10(30),log10(80)),size = 3)

# Figure 2C: within sample coincidence ####
# data
wsc <- read.csv("data/pc_withindonor_beta_down-sampled.csv") %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))

# stats
stats <- wsc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot 
p2C <- ggplot(wsc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "blue") +
  My_Theme+
  labs(x = "Sample", 
       y = "Within donor co-incidence") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top",
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats, label = "p.adj.signif", size = 3)+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0))

# Figure 2D: cross-sample coincidence ####
# data
csc <- read.csv("data/pc_crossdonor_beta_down-sampled.csv") %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))
# stats
stats <- csc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot
p2D <- ggplot(csc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "blue") +
  My_Theme+
  labs(x = "Sample", 
       y = "Between donor co-incidence") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top",
        axis.title.x = element_blank()) +
  stat_pvalue_manual(stats, label = "p.adj.signif", size = 3)+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0))

# Figure 2E: TST_D7 cross-sample coincidence stratified by HLA overlap ####
# HLA overlap
a <- read.csv("data/hladist_TST_D7_beta_down-sampled_mhcII.csv") %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class II")
b <- read.csv("data/hladist_TST_D7_beta_down-sampled_mhcI.csv") %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class I")
hla <- a %>% full_join(b)

# coincidence 
c <- read.csv("data/tcrsharingprob_TST_D7_beta_down-sampled_mhcII.csv") %>%
  rename("sharing.prob" = X0) %>%
  mutate(mhc = "MHC class II")
d <- read.csv("data/tcrsharingprob_TST_D7_beta_down-sampled_mhcI.csv") %>%
  rename("sharing.prob" = X0) %>%
  mutate(mhc = "MHC class I")

csc <- c %>% full_join(d)

# combine data
dat <- full_join(hla,csc)

# regression coefficient + 95% CI
reg1 <- dat %>%
  filter(mhc == "MHC class I") %>%
  lm(sharing.prob~HLA.dist, data=.)
coeff1 <- formatC(reg1$coefficients[2],format="e",digits=2)
ci1 <- formatC(confint(reg1, 'HLA.dist', level=0.95),format="e",digits=2)

reg2 <- dat %>%
  filter(mhc == "MHC class II") %>%
  lm(sharing.prob~HLA.dist, data=.)
coeff2 <- formatC(reg2$coefficients[2],format="e",digits=2)
ci2 <- formatC(confint(reg2, 'HLA.dist', level=0.95),format="e",digits=2)

dat_text <- data.frame(
  label = c(paste0("\u03B2","=",coeff1," (",ci1[1],"-",ci1[2],")"),
            paste0("\u03B2","=",coeff2," (",ci2[1],"-",ci2[2],")")),
  mhc = c('MHC class I','MHC class II'),
  x = c(2.5,4),
  y= c(1e-05,1e-05)
)

# plot
p2E <- ggplot(dat,aes(x=HLA.dist,y=sharing.prob)) +
  geom_jitter(aes(alpha = 0.1), colour="black", size=0.5) +
  geom_smooth(method = "lm")+
  scale_y_log10(limits = c(1e-8,1e-5))+
  facet_wrap(~mhc,scales = "free_x")+
  labs(x = "HLA overlap", 
       y = "Between donor co-incidence")+
  My_Theme+
  theme(legend.position="none") +
  geom_text(data = dat_text,
            mapping = aes(x=x,y=y,label=label),
            size = 6/.pt)

# assemble figure ####
p3 <- ggarrange(p2C,p2D,p2E,
                nrow = 1,
                widths = c(1,1,2),
                labels = list("C","D","E"),
                font.label = list(size = 10, face = "bold", colour = "black"))
ggarrange(p2A,p2B,p3,
          nrow=3,
          heights = c(1,1.3,1.1),
          labels = list("A","B"),
          font.label = list(size = 10, face = "bold", colour = "black"))
ggsave("Figure2.svg", 
       units = "cm", width = 17, height =22 , dpi=300)
