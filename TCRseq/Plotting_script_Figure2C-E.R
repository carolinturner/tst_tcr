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

# Figure 2C: within sample coincidence ####
wsc <- read.csv("data/pc_withindonor_beta_down-sampled.csv")
wsc <- wsc %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))

stats <- wsc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot (save as svg 300x350)
ggplot(wsc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "blue") +
  My_Theme+
  labs(x = "Sample", 
       y = "TCR sharing probability") +
  ggtitle("Within donor") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top")+
  stat_pvalue_manual(stats, label = "p.adj.signif")+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0))


# Figure 2D: cross-sample coincidence ####
csc <- read.csv("data/pc_crossdonor_beta_down-sampled.csv")

csc <- csc %>%
  mutate(tissue = recode(tissue,
                         'TST_D2'= "Day 2 TST",
                         'TST_D7'= "Day 7 TST"))
stats <- csc %>%
  pairwise_wilcox_test(data = ., logpc~tissue, p.adjust.method = "fdr") %>%
  add_xy_position() 

# plot (save as svg 300x350)
ggplot(csc, aes(x=tissue, y=logpc))+
  geom_boxplot(colour = "blue") +
  My_Theme+
  labs(x = "Sample", 
       y = "TCR sharing probability") +
  ggtitle("Cross-donor") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        legend.position = "right", legend.justification = "top") +
  stat_pvalue_manual(stats, label = "p.adj.signif")+
  scale_y_continuous(labels = \(x) formatC(10^x,format="e", digits =0))

# Figure 2E: TST_D7 cross-sample coincidence stratified by HLA overlap ####
# load HLA overlap
a <- read.csv("data/hladist_TST_D7_beta_down-sampled_mhcII.csv")
a <- a %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class II")
b <- read.csv("data/hladist_TST_D7_beta_down-sampled_mhcI.csv")
b <- b %>%
  rename("HLA.dist" = X0) %>%
  mutate(mhc = "MHC class I")
hla <- a %>% full_join(b)

# load coincidence data
c <- read.csv("data/tcrsharingprob_TST_D7_beta_down-sampled_mhcII.csv")
c <- c %>%
  rename("sharing.prob" = X0) %>%
  mutate(mhc = "MHC class II")
d <- read.csv("data/tcrsharingprob_TST_D7_beta_down-sampled_mhcI.csv")
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

# summary
reg <- as.data.frame(rbind(ci2,ci1))
rownames(reg) <- NULL
reg$mhc <- c("MHC class I","MHC class II")
reg$beta <- c(coeff2,coeff1)
reg$beta <- formatC(reg$beta, format = "e", digits = 2)
reg$`2.5 %`<- formatC(reg$`2.5 %`,format='e', digits=2)
reg$`97.5 %`<- formatC(reg$`97.5 %`,format = 'e', digits=2)
reg <- reg %>%
  mutate("95% CI" = paste0("(",`2.5 %`," - ",`97.5 %`,")")) %>%
  select(mhc,beta,'95% CI')
write.csv(reg,"data/Figure2E_regression-analysis.csv",row.names = F)

# plot (save as svg 600x500)
ggplot(dat,aes(x=HLA.dist,y=sharing.prob)) +
  geom_jitter(aes(alpha = 0.5)) +
  geom_smooth(method = "lm")+
  scale_y_log10(limits = c(1e-8,1e-5))+
  scale_x_continuous(breaks = seq(min(dat$HLA.dist), max(dat$HLA.dist), by = 1))+
  facet_wrap(~mhc,scales = "free_x")+
  labs(x = "HLA overlap", 
       y = "TCR sharing probability (cross-donor)")+
  My_Theme+
  theme(legend.position="none")
