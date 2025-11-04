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
  plot.margin = unit(c(0.5,0.25,0,0.2),"cm"),
  panel.spacing.x = unit(0.05,"cm")
)

# Figure S6A-C: within sample coincidence ####

## alpha, down-sampled
wsc <- read.csv("data/pc_withindonor_alpha_down-sampled.csv")
wsc <- wsc %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))

stats <- wsc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot
pS6A <- ggplot(wsc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "red") +
  My_Theme+
  labs(x = "Sample", 
       y = "Within donor co-incidence")+ 
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top",
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats, label = "p.adj.signif", size = 2)+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0),
                     expand = expansion(mult = c(0.05, 0.15)))

## beta, full repertoires
wsc <- read.csv("data/pc_withindonor_beta_full-repertoires.csv")
wsc <- wsc %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))

stats <- wsc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot
pS6B <- ggplot(wsc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "blue") +
  My_Theme+
  labs(x = "Sample", 
       y = "Within donor co-incidence")+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top",
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats, label = "p.adj.signif", size = 2)+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0),
                     expand = expansion(mult = c(0.05, 0.15)))


## alpha, full repertoires
wsc <- read.csv("data/pc_withindonor_alpha_full-repertoires.csv")
wsc <- wsc %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))

stats <- wsc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot
pS6C <- ggplot(wsc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "red") +
  My_Theme+
  labs(x = "Sample", 
       y = "Within donor co-incidence") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top",
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats, label = "p.adj.signif", size = 2)+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0),
                     expand = expansion(mult = c(0.05, 0.15)))


# Figure S6D-F: cross-sample coincidence ####

## alpha, down-sampled
csc <- read.csv("data/pc_crossdonor_alpha_down-sampled.csv")

csc <- csc %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))
stats <- csc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot
pS6D <- ggplot(csc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "red") +
  My_Theme+
  labs(x = "Sample", 
       y = "Between donor co-incidence") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top",
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats, label = "p.adj.signif", size = 2)+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0),
                     expand = expansion(mult = c(0.05, 0.15)))


## beta, full repertoires
csc <- read.csv("data/pc_crossdonor_beta_full-repertoires.csv")

csc <- csc %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))
stats <- csc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot
pS6E <- ggplot(csc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "blue") +
  My_Theme+
  labs(x = "Sample", 
       y = "Between donor co-incidence") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top",
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats, label = "p.adj.signif", size = 2)+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0),
                     expand = expansion(mult = c(0.05, 0.15)))


## alpha, full repertoires
csc <- read.csv("data/pc_crossdonor_alpha_full-repertoires.csv")

csc <- csc %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))
stats <- csc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot
pS6F <- ggplot(csc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "red") +
  My_Theme+
  labs(x = "Sample", 
       y = "Between donor co-incidence") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top",
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats, label = "p.adj.signif", size = 2)+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0),
                     expand = expansion(mult = c(0.05, 0.15)))


# Figure S6G-I: TST_D7 cross-sample coincidence stratified by HLA overlap ####

## alpha, down-sampled
# HLA overlap
a <- read.csv("data/hladist_TST_D7_alpha_down-sampled_mhcII.csv")
a <- a %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class II")
b <- read.csv("data/hladist_TST_D7_alpha_down-sampled_mhcI.csv")
b <- b %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class I")
hla <- a %>% full_join(b)

# coincidence
c <- read.csv("data/tcrsharingprob_TST_D7_alpha_down-sampled_mhcII.csv")
c <- c %>%
  rename("sharing.prob" = X0) %>%
  mutate(mhc = "MHC class II")
d <- read.csv("data/tcrsharingprob_TST_D7_alpha_down-sampled_mhcI.csv")
d <- d %>%
  rename("sharing.prob" = X0) %>%
  mutate(mhc = "MHC class I")

csc <- c %>% full_join(d)

# combine data
dat <- full_join(hla,csc)

# regression coefficient + 95% CI
reg1 <- dat %>%
  filter(mhc == "MHC class I") %>%
  lm(sharing.prob~HLA.dist, data=.)
coeff1 <- reg1$coefficients[2]
ci1 <- confint(reg1, 'HLA.dist', level=0.95)

reg2 <- dat %>%
  filter(mhc == "MHC class II") %>%
  lm(sharing.prob~HLA.dist, data=.)
coeff2 <- reg2$coefficients[2]
ci2 <- confint(reg2, 'HLA.dist', level=0.95)

# plot
pS6G <- ggplot(dat,aes(x=HLA.dist,y=sharing.prob)) +
  geom_jitter(aes(alpha = 0.5)) +
  geom_smooth(method = "lm")+
  scale_y_log10(limits = c(7e-7,3e-5))+
  facet_wrap(~mhc,scales = "free_x")+
  labs(x = "HLA overlap", 
       y = "Between donor co-incidence")+
  My_Theme+
  theme(legend.position="none")


## beta, full repertoires
# HLA overlap
a <- read.csv("data/hladist_TST_D7_beta_full-repertoires_mhcII.csv")
a <- a %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class II")
b <- read.csv("data/hladist_TST_D7_beta_full-repertoires_mhcI.csv")
b <- b %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class I")
hla <- a %>% full_join(b)

# coincidence
c <- read.csv("data/tcrsharingprob_TST_D7_beta_full-repertoires_mhcII.csv")
c <- c %>%
  rename("sharing.prob" = X0) %>%
  mutate(mhc = "MHC class II")
d <- read.csv("data/tcrsharingprob_TST_D7_beta_full-repertoires_mhcI.csv")
d <- d %>%
  rename("sharing.prob" = X0) %>%
  mutate(mhc = "MHC class I")

csc <- c %>% full_join(d)

# combine data
dat <- full_join(hla,csc)

# regression coefficient + 95% CI
reg1 <- dat %>%
  filter(mhc == "MHC class II") %>%
  lm(sharing.prob~HLA.dist, data=.)
coeff1 <- reg1$coefficients[2]
ci1 <- confint(reg1, 'HLA.dist', level=0.95)

reg2 <- dat %>%
  filter(mhc == "MHC class I") %>%
  lm(sharing.prob~HLA.dist, data=.)
coeff2 <- reg2$coefficients[2]
ci2 <- confint(reg2, 'HLA.dist', level=0.95)

# plot
pS6H <- ggplot(dat,aes(x=HLA.dist,y=sharing.prob)) +
  geom_jitter(aes(alpha = 0.5)) +
  geom_smooth(method = "lm")+
  scale_y_log10(limits = c(1e-8,2e-5))+
  facet_wrap(~mhc,scales = "free_x")+
  labs(x = "HLA overlap", 
       y = "Between donor co-incidence")+
  My_Theme+
  theme(legend.position="none")


## alpha, full repertoires
# HLA overlap
a <- read.csv("data/hladist_TST_D7_alpha_full-repertoires_mhcII.csv")
a <- a %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class II")
b <- read.csv("data/hladist_TST_D7_alpha_full-repertoires_mhcI.csv")
b <- b %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class I")
hla <- a %>% full_join(b)

# coincidence
c <- read.csv("data/tcrsharingprob_TST_D7_alpha_full-repertoires_mhcII.csv")
c <- c %>%
  rename("sharing.prob" = X0) %>%
  mutate(mhc = "MHC class II")
d <- read.csv("data/tcrsharingprob_TST_D7_alpha_full-repertoires_mhcI.csv")
d <- d %>%
  rename("sharing.prob" = X0) %>%
  mutate(mhc = "MHC class I")

csc <- c %>% full_join(d)

# combine data
dat <- full_join(hla,csc)

# regression coefficient + 95% CI
reg1 <- dat %>%
  filter(mhc == "MHC class II") %>%
  lm(sharing.prob~HLA.dist, data=.)
coeff1 <- reg1$coefficients[2]
ci1 <- confint(reg1, 'HLA.dist', level=0.95)

reg2 <- dat %>%
  filter(mhc == "MHC class I") %>%
  lm(sharing.prob~HLA.dist, data=.)
coeff2 <- reg2$coefficients[2]
ci2 <- confint(reg2, 'HLA.dist', level=0.95)

# plot
pS6I <- ggplot(dat,aes(x=HLA.dist,y=sharing.prob)) +
  geom_jitter(aes(alpha = 0.5)) +
  geom_smooth(method = "lm")+
  scale_y_log10(limits = c(2e-7,2e-4))+
  facet_wrap(~mhc,scales = "free_x")+
  labs(x = "HLA overlap", 
       y = "Between donor co-incidence")+
  My_Theme+
  theme(legend.position="none")

# assemble figure ####
ggarrange(pS6A,pS6B,pS6C,pS6D,pS6E,pS6F,pS6G,pS6H,pS6I,
          nrow=3,
          ncol=3,
          labels = list("A","B","C","D","E","F","G","H","I"),
          font.label = list(size = 10, face = "bold", colour = "black"))
ggsave("figures/FigureS6.svg", 
       units = "cm", width = 17, height =17 , dpi=300)

