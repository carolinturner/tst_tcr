library(tidyverse)
library(ggh4x)

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

# Figure S7A ####
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

# Figure S7B ####
# remove clusters with odds ratio = 0 as can't log10 transform
# replace odds ratio = Inf with very big, numerical value to allow plotting
d1 <- read.csv("output/metaclonotypist_alpa_mhcII/clusterassociation_alpha_mhcII.csv") %>%
  filter(odds_ratio > 0) %>%
  mutate(odds_ratio = replace(odds_ratio, odds_ratio == "Inf", 400),
         log_p = -log10(pvalue),
         dataset = "True data",
         mhc = "MHC-II") %>%
  select(odds_ratio,log_p,dataset,mhc,significant)
d2 <- read.csv("output/metaclonotypist_alpha_mhcI/clusterassociation_alpha_mhcI.csv") %>%
  filter(odds_ratio > 0) %>%
  mutate(odds_ratio = replace(odds_ratio, odds_ratio == "Inf", 400),
         log_p = -log10(pvalue),
         dataset = "True data",
         mhc = "MHC-I") %>%
  select(odds_ratio,log_p,dataset,mhc,significant)
d3 <- read.csv("output/metaclonotypist_alpha_mhcII/clusterassociation_shuffled_alpha_mhcII.csv") %>%
  filter(odds_ratio > 0) %>%
  mutate(odds_ratio = replace(odds_ratio, odds_ratio == "Inf", 400),
         log_p = -log10(pvalue),
         dataset = "Shuffled data",
         mhc = "MHC-II") %>%
  select(odds_ratio,log_p,dataset,mhc,significant)
d4 <- read.csv("output/metaclonotypist_alpha_mhcI/clusterassociation_shuffled_alpha_mhcI.csv") %>%
  filter(odds_ratio > 0) %>%
  mutate(odds_ratio = replace(odds_ratio, odds_ratio == "Inf", 400),
         log_p = -log10(pvalue),
         dataset = "Shuffled data",
         mhc = "MHC-I") %>%
  select(odds_ratio,log_p,dataset,mhc,significant)

dat <- rbind(d1,d2,d3,d4)
dat$dataset <- factor(dat$dataset,levels=c("True data","Shuffled data"))
dat$mhc <- factor(dat$mhc, levels=c("MHC-II","MHC-I"))

p <- ggplot(dat, aes(x=odds_ratio, y=log_p, colour=significant))+
  geom_point(alpha=0.5)+
  facet_grid(dataset~mhc)+
  scale_colour_manual(values = c("darkgrey","darkblue"))+
  scale_x_log10()+
  My_Theme +
  labs(y= "-log10 p-value",
       x = "odds ratio")
ggsave("FigureS7B.png",p, units = "cm", width = 15.9, height = 11.9, dpi = 600)

# Figure S7C ####
# load metaclonotypist output
b1 <- read.csv("output/metaclonotypist_alpha_mhcI/hlametaclonotypes_alpha_mhcI.csv")
b2 <- read.csv("output/metaclonotypist_alpha_mhcII/hlametaclonotypes_alpha_mhcII.csv")
# alternative (using Supplementary tables)
#b1 <- read.csv("data/TableS8.csv")
#b2 <- read.csv("data/TableS5.csv")

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

# plot (save as svg 900x400)
ggplot(b,aes(x=reorder(hla,desc(hla)),y=n)) +
  geom_bar(stat = "identity",fill="red") +
  facet_wrap(~dataset,scale="free_x") +
  force_panelsizes(cols = c(15,1)) +
  My_Theme +
  labs(x="HLA association",
       y="Number of metaclones") +
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0))
