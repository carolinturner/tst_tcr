library(tidyverse)
library(scales)
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
  panel.spacing.x = unit(0.05,"cm")
)

# Figure 5A ####
dat <- read.csv("data/Figure5A_odds-ratio.csv") %>%
  mutate(Log10.OR=log10(oddsratio),
         Log10.ORLL=log10(CI_lower),
         Log10.ORUL=log10(CI_higher),
         DatasetLabel = recode(dataset,
                               "PBMC" = "PPD:TT stimulated.PBMC",
                               "Tcells" = "Mtb:SARS-CoV2 reactive.Tcells",
                               "scLung" = "TB:Cancer Lung",
                               "Blood" = "TB:Cancer Blood",
                               "CD4-T" = "Lung:Blood CD4.Tcells",
                               "Lung" = "TB:Cancer Lung"),
         Seq.Method = recode(dataset,
                             "PBMC" = "Bulk",
                             "Tcells" = "SingleCell",
                             "scLung" = "SingleCell",
                             "Blood" = "Bulk",
                             "CD4-T" = "Bulk",
                             "Lung" = "Bulk"))
dat$algorithm <- factor(dat$algorithm , levels = c("Discovery TST CDR3s",
                                                   "Metaclonotypist regex",
                                                   "Gliph2 pattern",
                                                   "Gliph2-matching CDR3s"))
dat$DatasetLabel <- factor(dat$DatasetLabel, levels = c("PPD:TT stimulated.PBMC",
                                                        "Mtb:SARS-CoV2 reactive.Tcells",
                                                        "TB:Cancer Lung",
                                                        "TB:Cancer Blood",
                                                        "Lung:Blood CD4.Tcells"))

p5A <- ggplot(dat, aes(x=DatasetLabel,y=Log10.OR,colour = algorithm, shape=Seq.Method)) +
  geom_point(stat = "identity", position = position_dodge2(width=0.5), size=2, alpha=0.8)+
  geom_errorbar(aes(ymin = Log10.ORLL,ymax = Log10.ORUL), width=0.5, position=position_dodge2(width=0.5))+
  scale_x_discrete(labels = label_wrap(10))+
  scale_colour_manual(values=c("darkgrey", "#0F2080","#ac6914ff", "#ac141bff"))+
  labs(x="Dataset",
       y="Log10 Odds Ratio",
       colour="TCR set",
       shape="Sequencing method")+
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        plot.margin = unit(c(0.2,0,0.5,0.5),"cm"))+
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2))

# Figure 5B ####
a <- read.csv("data/summary_metaclone-abundance_down-sampled_beta_expanded_gr0.csv")
b <- read.csv("data/summary_metaclone-abundance_down-sampled_beta_expanded_gr1.csv")

# beta chain data
c <- a %>% mutate(Clone.Size = "All TCRs")
d <- b %>% mutate(Clone.Size = "Expanded TCRs")

summary <- rbind(c,d) %>%
  mutate(tissue = recode(tissue,
                         TST_D2 = "Day 2 TST",
                         TST_D7 = "Day 7 TST"))
# stats
stats.all <- summary %>%
  group_by(Clone.Size) %>%
  pairwise_wilcox_test(data = .,mc.pct~tissue, p.adjust.method = "fdr")

# plot
p5B <- ggplot(summary, aes(x=tissue,y=mc.pct))+
  geom_boxplot(colour="blue")+
  facet_wrap(~Clone.Size,ncol = 5)+
  scale_y_log10(limits = c(0.03,100))+
  labs(y="% of TCRs (down-sampled)") +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.2,0.5,0.5,0.5),"cm"))+
  stat_pvalue_manual(stats.all, label = "p.adj.signif",y.position=c(log10(15),log10(30),log10(60)), size = 3)

# Figure 5C ####
# in vitro vs mc
keep <- c("Sample","UIN","pct.PPD.only","pct.mc")

# down-sampled beta: Day 2 TST
a <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr0.csv") %>% select(all_of(keep))
b <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr1.csv") %>% select(all_of(keep))
c <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr2.csv") %>% select(all_of(keep))
d <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr3.csv") %>% select(all_of(keep))
e <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr4.csv") %>% select(all_of(keep))

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e)

beta.d2 <- summary %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public")) %>% 
  mutate(Sample="Day 2 TST")

# down-sampled beta: Day 7 TST
a <- read.csv("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr0.csv") %>% select(all_of(keep))
b <- read.csv("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr1.csv") %>% select(all_of(keep))
c <- read.csv("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr2.csv") %>% select(all_of(keep))
d <- read.csv("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr3.csv") %>% select(all_of(keep))
e <- read.csv("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr4.csv") %>% select(all_of(keep))

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e)

beta.d7 <- summary %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public"))%>% 
  mutate(Sample="Day 7 TST")

# combine D2 and D7 data
beta.all <- rbind(beta.d2,beta.d7) %>% 
  mutate(TST.day= as.character(ifelse(Sample=="Day 2 TST", 2, 7)))

p5C <- ggplot(beta.all, aes(x=TST.day, y=percentage, colour=TCR.class))+
  facet_wrap(~Clone.Size, nrow=1)+
  geom_point(size=1, alpha=0.5) +
  geom_line(aes(group=interaction(UIN,TCR.class), 
                colour = TCR.class), alpha=0.5)+
  scale_colour_manual(values = c("orange","navy"))+
  scale_y_continuous(transform = "log10")+
  labs(y="% of TCRs (down-sampled)",
       x = "TST sample day",
       colour="TCR class",
       title="TCR clone size")+
  My_Theme +
  theme(plot.margin = unit(c(0.2,0,0.5,0.5),"cm"))

# Figure 5D ####
a <- read.csv("data/Publicity_all-vs-published-vs-metaclone_down-sampled_beta_expanded_gr0.csv")
b <- read.csv("data/Publicity_all-vs-published-vs-metaclone_down-sampled_beta_expanded_gr1.csv")

a <- a %>% mutate(Clone.Size = "All TCRs")
b <- b %>% mutate(Clone.Size = "Expanded TCRs")

summary <- rbind(a,b) %>%
  dplyr::rename(TCR = "dataset") %>%
  mutate(TCR = recode(TCR,
                      "CDR3 with published Mtb reactivity" = "Published Mtb CDR3s",
                      "CDR3" = "D7 TST CDR3s",
                      "metaclone"="D7 TST metaclones"))

# plot
p5D <- ggplot(summary, aes(x=rank,y=cum.prop.people,color=TCR))+
  geom_line() +
  geom_point(alpha=0.8) +
  scale_colour_manual(values=c("darkred","navy","orange")) +
  facet_wrap(~Clone.Size)+
  scale_x_log10()+
  labs(x="Number of TCRs (ranked by publicity)",
       y="Proportion of participants \n(cumulative)") +
  My_Theme +
  theme(legend.position = "right",legend.justification="top",
        panel.spacing = unit(0.5, "cm", data = NULL),
        plot.margin = unit(c(0.2,1,0.5,0.5),"cm"))

# Figure 5E ####
mc <- read.csv("data/TableS4.csv") %>%
  mutate(hla.pct = count_allele/(count_allele+count_other)*100,
         hla.freq = count_allele+count_other)

p5E <- mc %>%
  mutate_at("index", as.character) %>% 
  ggplot()+
  geom_point(aes(x=hla.freq,y=hla.pct,colour=index))+
  coord_cartesian(ylim = c(0,100)) +
  labs(x="HLA frequency",
       y="Percentage of all participants \nwith cognate HLA",
       title = "Individual metaclones by colour")+
  My_Theme +
  theme(plot.margin = unit(c(0.2,0.5,0.5,0),"cm"),
        legend.position = "none")
p5E

# assemble figure ####
r2 <- ggarrange(p5B,p5C,
                ncol = 2,
                widths = c(1,1.2),
                labels = c("B","C"),
                font.label = list(size = 10, face = "bold", colour = "black"))
r3 <- ggarrange(p5D,p5E,
                ncol = 2,
                widths = c(2,1),
                labels = c("D","E"),
                font.label = list(size = 10, face = "bold", colour = "black"))

ggarrange(p5A,r2,r3,
          nrow = 3,
          labels = c("A","",""),
          font.label = list(size = 10, face = "bold", colour = "black"))

ggsave("Figure5.svg", 
       units = "cm", width = 17, height =20 , dpi=300)
