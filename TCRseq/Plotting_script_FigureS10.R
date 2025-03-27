library(tidyverse)
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
  panel.spacing.x = unit(0.05,"cm"),
  plot.margin = unit(c(0.5,0.6,0,0.2),"cm")
)

# Figure S10A: private vs public TST enrichment (full repertoires, beta) ####
keep <- c("Sample","UIN","pct.PPD.only","pct.mc")

## Day 2 TST 
a <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr0.csv") %>% select(all_of(keep))
b <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr1.csv") %>% select(all_of(keep))
c <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr2.csv") %>% select(all_of(keep))
d <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr3.csv") %>% select(all_of(keep))
e <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr4.csv") %>% select(all_of(keep))

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e)

d2 <- summary %>%
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public")) %>% 
  mutate(Sample="Day 2 TST")

## Day 7 TST 
a <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr0.csv") %>% select(all_of(keep))
b <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr1.csv") %>% select(all_of(keep))
c <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr2.csv") %>% select(all_of(keep))
d <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr3.csv") %>% select(all_of(keep))
e <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr4.csv") %>% select(all_of(keep))

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e)

d7 <- summary %>%
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public")) %>% 
  mutate(Sample="Day 7 TST")

## combine D2 and D7 data
all <- rbind(d2,d7) %>% 
  mutate(TST.day= as.character(ifelse(Sample=="Day 2 TST", 2, 7)))

# plot
pS10A <- ggplot(all, aes(x=TST.day, y=percentage, colour=TCR.class))+
  facet_wrap(~Clone.Size, nrow=1)+
  geom_point(size=1, alpha=0.5) +
  geom_line(aes(group=interaction(UIN,TCR.class), 
                colour = TCR.class), alpha=0.5)+
  scale_colour_manual(values = c("orange","navy"))+
  scale_y_continuous(transform = "log10")+
  labs(y="% of TCRs (full repertoire)",
       x = "TST sample day",
       colour="TCR class",
       title="TCR clone size")+
  My_Theme

# Figure S10B: expansion d2-d7 private vs public (down-sampled, beta) ####
keep <- c("Sample","UIN","pct.PPD.only","pct.mc")

## down-sampled beta: Day 2 TST 
a <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr0.csv")
b <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr1.csv")
c <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr2.csv")
d <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr3.csv")
e <- read.csv("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr4.csv")

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e)

# replace mc count of 0 with count of 1 and recalculate pct
summary_new <- summary %>%
  mutate(number.mc = ifelse(number.mc == 0, 1, number.mc),
         pct.mc = number.mc/total.TCRs*100) %>%
  select(all_of(keep),Clone.Size)

d2 <- summary_new %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public"))

## down-sampled beta: Day 7 TST
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

d7 <- summary %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public"))

## delta analysis 
d2 <- d2 %>% 
  select(UIN, Clone.Size, TCR.class, percentage) %>% 
  rename(d2.percentage=percentage)
d7 <- d7 %>% 
  select(UIN, Clone.Size, TCR.class, percentage)%>% 
  rename(d7.percentage=percentage)
delta <- d2 %>% 
  left_join(d7, by=c("UIN", "Clone.Size", "TCR.class")) %>% 
  mutate(FC.percentage=d7.percentage/d2.percentage)

# stats
stats <- delta %>% 
  group_by(Clone.Size) %>% 
  wilcox_test(FC.percentage~TCR.class) %>% 
  add_xy_position()

# plot
pS10B <- ggplot(delta, aes(x=TCR.class,y=FC.percentage,colour=TCR.class))+
  geom_boxplot(outliers=F)+
  facet_wrap(~Clone.Size, nrow=1)+
  scale_colour_manual(values = c("orange","navy"))+
  stat_pvalue_manual(stats,label="p={p}", y.position = 20, size=3, tip.length=0.002)+
  scale_y_continuous(limits=c(0,21))+
  labs(x="PPD-reactive TCR class",
       y="Fold change in % TCRs \n(down-sampled))",
       title="TCR clone size")+
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -30, hjust = 0),
        legend.position="none")

# Figure S10C: expansion d2-d7 private vs public (full repertoires, beta) ####
keep <- c("Sample","UIN","pct.PPD.only","pct.mc")

## down-sampled beta: Day 2 TST 
a <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr0.csv")
b <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr1.csv")
c <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr2.csv")
d <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr3.csv")
e <- read.csv("data/invitro-vs-mc_results_D2_full-repertoires_expanded_gr4.csv")

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e)

# replace mc count of 0 with count of 1 and recalculate pct
summary_new <- summary %>%
  mutate(number.mc = ifelse(number.mc == 0, 1, number.mc),
         pct.mc = number.mc/total.TCRs*100) %>%
  select(all_of(keep),Clone.Size)

d2 <- summary_new %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public"))

## down-sampled beta: Day 7 TST
a <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr0.csv") %>% select(all_of(keep))
b <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr1.csv") %>% select(all_of(keep))
c <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr2.csv") %>% select(all_of(keep))
d <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr3.csv") %>% select(all_of(keep))
e <- read.csv("data/invitro-vs-mc_results_D7_full-repertoires_expanded_gr4.csv") %>% select(all_of(keep))

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e)

d7 <- summary %>% 
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public"))

# delta analysis ####
d2 <- d2 %>% 
  select(UIN, Clone.Size, TCR.class, percentage) %>% 
  rename(d2.percentage=percentage)
d7 <- d7 %>% 
  select(UIN, Clone.Size, TCR.class, percentage)%>% 
  rename(d7.percentage=percentage)
delta <- d2 %>% 
  left_join(d7, by=c("UIN", "Clone.Size", "TCR.class")) %>% 
  mutate(FC.percentage=d7.percentage/d2.percentage)

# stats
stats <- delta %>% 
  group_by(Clone.Size) %>% 
  wilcox_test(FC.percentage~TCR.class) %>% 
  add_xy_position()

# plot
pS10C <- ggplot(delta, aes(x=TCR.class,y=FC.percentage, colour = TCR.class))+
  geom_boxplot(outliers=F)+
  facet_wrap(~Clone.Size, nrow=1)+
  scale_colour_manual(values = c("orange","navy"))+
  stat_pvalue_manual(stats,label="p={p}", y.position = 20, size=3, tip.length=0.01)+
  scale_y_continuous(limits=c(0,21))+
  labs(x="PPD-reactive TCR class",
       y="Fold change in % TCRs \n(full repertoire))",
       title="TCR clone size")+
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -30, hjust = 0),
        legend.position="none")

# assemble figure ####
r1 <- ggarrange(pS10A,
                ncol = 2,
                widths = c(1.4,1),
                labels = c("A"),
                font.label = list(size = 10, face = "bold", colour = "black"))
r2 <- ggarrange(pS10B,pS10C,
                labels = c("B","C"),
                font.label = list(size = 10, face = "bold", colour = "black"))
ggarrange(r1,r2,
          nrow = 2)
ggsave("FigureS10.svg", 
       units = "cm", width = 17, height =14 , dpi=300)
