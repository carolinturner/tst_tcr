library(tidyverse)

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

dat <- read.csv("data/Figure5A_percentage.csv")

# set levels
dat$algorithm <- factor(dat$algorithm,
                       levels = c("Metaclonotypist","Gliph2"))
dat$dataset <- factor(dat$dataset,
                     levels = c("PBMC (bulk-TCRseq)",
                                "T-cells (sc-TCRseq)",
                                "Lung (sc-TCRseq)",
                                "Lung (bulk-TCRseq)",
                                "Blood (bulk-TCRseq)",
                                "CD4-T (bulk-TCRseq)"))
dat$Stimulant <- factor(dat$Stimulant,
                       levels = c("PPD","TT",
                                  "Mtb","SARS-CoV2",
                                  "TB","Cancer",
                                  "TB lung","TB blood"))

# plot (save as svg 1150x550, re-format in Inkscape)
ggplot(dat, aes(x=Stimulant,y=mc.pct.total))+
  geom_col()+
  facet_wrap(algorithm ~ dataset, scales = "free",ncol=6)+
  My_Theme+
  labs(x="Stimulant/Disease/Tissue",
       y="% of dataset identified as match")
