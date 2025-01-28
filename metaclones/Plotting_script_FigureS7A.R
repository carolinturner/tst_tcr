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
  panel.background = element_rect(fill = "gray97"),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  legend.position = "right", legend.justification = "top"
)

# load data
dat <- read.csv("data/parameter_sweep.csv")
b1 <- dat %>%
  filter(chain == "alpha" & mincount == 2 & clustering == "leiden") %>%
  select(max_tcrdist,nmetaclones,sig_clonotype_fraction,id_fraction) %>%
  mutate(dataset = "true")
b2 <- dat %>%
  filter(chain == "alpha" & mincount == 2 & clustering == "leiden") %>%
  select(max_tcrdist,nmetaclones_shuffled,sig_clonotype_fraction_shuffled,id_fraction_shuffled) %>%
  mutate(dataset = "shuffled") %>%
  rename(nmetaclones = nmetaclones_shuffled,
         sig_clonotype_fraction = sig_clonotype_fraction_shuffled,
         id_fraction = id_fraction_shuffled)
b <- rbind(b1,b2)

# prep for plotting
b.long <- b %>%
  pivot_longer(cols = nmetaclones:id_fraction,
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         nmetaclones = "# HLA-metaclones",
                         sig_read_fraction = "% TCRs",
                         id_fraction = "% Participants"))

b.long$dataset <- factor(b.long$dataset, levels = c(
  "true",
  "shuffled"
))

# plot (save as svg 610x250)
p <- ggplot(b.long,aes(max_tcrdist,Value,colour=dataset))+
  geom_point()+
  geom_line()+
  scale_colour_manual(values=c("blue","orange"))+
  facet_wrap(~Metric,ncol=3,scales = "free_y")+
  labs(x = "TCRdist threshold",
       y = "")+
  My_Theme 