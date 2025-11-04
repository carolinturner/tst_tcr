library(tidyverse)
library(scales)
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

dat <- read.csv("data/Figure10B_odds-ratio_gliph-filters.csv") %>%
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
dat$algorithm <- factor(dat$algorithm , levels = c("Gliph2 + Filter 1",
                                                   "Gliph2 + Filter 1 + Filter 2 (p<0.05)",
                                                   "Gliph2 + Filter 1 + Filter 2 (FDR<0.1)",
                                                   "Gliph2 + Filter 2 (FDR<0.1)"))
dat$DatasetLabel <- factor(dat$DatasetLabel, levels = c("PPD:TT stimulated.PBMC",
                                                        "Mtb:SARS-CoV2 reactive.Tcells",
                                                        "TB:Cancer Lung",
                                                        "TB:Cancer Blood",
                                                        "Lung:Blood CD4.Tcells"))

p <- ggplot(dat, aes(x=DatasetLabel,y=Log10.OR,colour = algorithm, shape=Seq.Method)) +
  geom_point(stat = "identity", position = position_dodge2(width=0.5), size=2, alpha=0.8)+
  geom_errorbar(aes(ymin = Log10.ORLL,ymax = Log10.ORUL), width=0.5, position=position_dodge2(width=0.5))+
  scale_x_discrete(labels = label_wrap(10))+
  scale_colour_manual(values=c("darkgrey","#0F2080","#ac6914ff", "#ac141bff"))+
  labs(x="Dataset",
       y="Log10 Odds Ratio",
       colour="TCR set",
       shape="Sequencing method")+
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0),
        plot.margin = unit(c(0.2,0,0.5,0.5),"cm"))+
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2))
p

ggarrange(p,
          labels = c("B"),
          font.label = list(size = 10, face = "bold", colour = "black"))
ggsave("figures/FigureS10B.svg", 
       units = "cm", width = 17, height =7 , dpi=300)
