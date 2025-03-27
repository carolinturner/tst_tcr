library(tidyverse)
library(ggh4x)
library(ggpubr)

#My_Theme
t = 8 #size of text
m = 2 #size of margin around text
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
  legend.position = "top", legend.justification = "left",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0,0,0,0),
  legend.box.spacing = unit(c(0,0,0,0),"cm"),
  plot.margin = unit(c(0,0.2,0,0),"cm")
)

# Figure 4C ####
dat <- read.csv("data/parameter_sweep.csv")
b1 <- dat %>%
  filter(chain == "beta" & mincount == 2 & clustering == "leiden") %>%
  select(max_tcrdist,nmetaclones,sig_clonotype_fraction,id_fraction) %>%
  mutate(dataset = "true")
b2 <- dat %>%
  filter(chain == "beta" & mincount == 2 & clustering == "leiden") %>%
  select(max_tcrdist,nmetaclones_shuffled,sig_clonotype_fraction_shuffled,id_fraction_shuffled) %>%
  mutate(dataset = "shuffled") %>%
  rename(nmetaclones = nmetaclones_shuffled,
         sig_clonotype_fraction = sig_clonotype_fraction_shuffled,
         id_fraction = id_fraction_shuffled)
b <- rbind(b1,b2)
b <- b %>%
  mutate(sig_clonotype_fraction = 100*sig_clonotype_fraction,
         id_fraction = 100*id_fraction)

# prep for plotting
b.long <- b %>%
  pivot_longer(cols = nmetaclones:id_fraction,
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         nmetaclones = "Metaclones (N)",
                         sig_clonotype_fraction = "TCRs (%)",
                         id_fraction = "Participants (%)"))

b.long$dataset <- factor(b.long$dataset, levels = c(
  "true",
  "shuffled"
))

# plot
p4C <- ggplot(b.long,aes(max_tcrdist,Value,colour=dataset))+
  geom_point()+
  geom_line()+
  scale_colour_manual(values=c("blue","orange"))+
  facet_wrap(~Metric,ncol=3,scales = "free_y")+
  labs(x = "TCRdist threshold",
       y = "",
       col = "dataset")+
  My_Theme 
ggarrange(p4C,
          labels = list("C"),
          font.label = list(size = 10, face = "bold", colour = "black"))
ggsave("Figure4C.svg", 
       units = "cm", width = 12, height =6, dpi=300)

# Figure 4E ####
# load metaclonotypist output
b1 <- read.csv("output/metaclonotypist_beta_mhcI/hlametaclonotypes_beta_mhcI.csv")
b2 <- read.csv("output/metaclonotypist_beta_mhcII/hlametaclonotypes_beta_mhcII.csv")
# alternative (using Supplementary tables)
#b1 <- read.csv("data/TableS5.csv")
#b2 <- read.csv("data/TableS4.csv")


# calculate number of HLA associations
b1 <- b1 %>%
  group_by(hla) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  mutate(dataset = "MHC-I")
b1$hla <- factor(b1$hla, levels=b1$hla[order(b1$n)])

b2 <- b2 %>%
  group_by(hla) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  mutate(dataset = "MHC-II")
b2$hla <- factor(b2$hla, levels=b2$hla[order(b2$n)])

b <- rbind(b1,b2)
b$dataset <- factor(b$dataset,levels = c("MHC-II","MHC-I"))

# plot 
p4E <- ggplot(b,aes(x=reorder(hla,desc(hla)),y=n)) +
  geom_bar(stat = "identity",fill="blue") +
  facet_wrap(~dataset,scale="free_x") +
  force_panelsizes(cols = c(15,1)) +
  My_Theme +
  labs(x="HLA association",
       y="Metaclones (N)") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))
ggarrange(p4E,
          labels = list("E"),
          font.label = list(size = 10, face = "bold", colour = "black"))
ggsave("Figure4E.svg", 
       units = "cm", width = 17, height =6, dpi=300)
