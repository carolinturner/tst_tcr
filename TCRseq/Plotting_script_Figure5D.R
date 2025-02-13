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

a <- read.csv("data/Publicity_all-vs-published-vs-metaclone_down-sampled_beta_expanded_gr0.csv")
b <- read.csv("data/Publicity_all-vs-published-vs-metaclone_down-sampled_beta_expanded_gr1.csv")

a <- a %>% mutate(Clone.Size = "All CDR3s")
b <- b %>% mutate(Clone.Size = "Expanded CDR3s")

summary <- rbind(a,b) %>% rename(TCR = "dataset") %>%
  mutate(TCR = recode(TCR,
                      "CDR3 with published Mtb reactivity" = "CDR3 (published Mtb reactivity)",
                      "CDR3" = "CDR3 (D7 TST)"))

# plot (save as svg 600x400)
ggplot(summary, aes(x=rank,y=cum.prop.people,color=TCR))+
  geom_line() +
  geom_point() +
  scale_colour_manual(values=c("orange","darkred","navy")) +
  facet_wrap(~Clone.Size)+
  scale_x_log10()+
  labs(x="number of TCRs (ranked by publicity)",
       y="proportion of participants (cumulative)") +
  My_Theme +
  theme(legend.position = "top",legend.justification="left")
