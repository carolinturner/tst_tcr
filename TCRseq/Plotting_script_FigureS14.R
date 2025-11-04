library(tidyverse)
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


# full repertoires beta ####
a <- read.csv("data/Publicity_all-vs-published-vs-metaclone_full-repertoires_beta_expanded_gr0.csv")
b <- read.csv("data/Publicity_all-vs-published-vs-metaclone_full-repertoires_beta_expanded_gr1.csv")

a <- a %>% mutate(Clone.Size = "All TCRs")
b <- b %>% mutate(Clone.Size = "Expanded TCRs")

summary <- rbind(a,b) %>% rename(TCR = "dataset") %>%
  mutate(TCR = recode(TCR,
                      "CDR3 with published Mtb reactivity" = "CDR3 (published Mtb reactivity)",
                      "CDR3" = "CDR3 (D7 TST)"))

pS14 <- ggplot(summary, aes(x=rank,y=cum.prop.people,color=TCR))+
  geom_line() +
  geom_point() +
  scale_colour_manual(values=c("orange","darkred","navy")) +
  facet_wrap(~Clone.Size)+
  scale_x_log10()+
  labs(x="number of TCRs (ranked by publicity)",
       y="proportion of participants \n(cumulative)") +
  My_Theme +
  #theme(legend.position = "top",legend.justification="left") +
  ggtitle("beta (full repertoire)")

ggsave("figures/FigureS14.svg", 
       units = "cm", width = 12, height =6 , dpi=300)
