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

# in vitro vs mc
keep <- c("chain","Sample","UIN","total.CDR3s","number.PPD.only","pct.PPD.only","number.mc","pct.mc")

# down-sampled beta: Day 2 TST
a <- read.csv(paste0("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr0.csv")) %>%
  select(all_of(keep))
b <- read.csv(paste0("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr1.csv")) %>%
  select(all_of(keep))
c <- read.csv(paste0("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr2.csv")) %>%
  select(all_of(keep))
d <- read.csv(paste0("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr3.csv")) %>%
  select(all_of(keep))
e <- read.csv(paste0("data/invitro-vs-mc_results_D2_down-sampled_expanded_gr4.csv")) %>%
  select(all_of(keep))

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e)

beta <- summary %>% filter(chain == "beta") %>%
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public"))

# plot (save as svg 600x350)
ggplot(beta, aes(x=TCR.class,y=percentage,colour=TCR.class)) +
  geom_line(aes(group=UIN),color="black")+
  geom_point(size=2,position = position_dodge(width=0.7)) +
  facet_wrap(~Clone.Size,scales="free_x",strip.position = "bottom",ncol=5)+
  My_Theme +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0,units="line"))+
  scale_colour_manual(values = c("orange","navy")) +
  labs(y="% of total CDR3s (down-sampled)",
       x = "Clone size") +
  ggtitle("Day 2 TST")

# down-sampled beta: Day 7 TST
a <- read.csv(paste0("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr0.csv")) %>%
  select(all_of(keep))
b <- read.csv(paste0("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr1.csv")) %>%
  select(all_of(keep))
c <- read.csv(paste0("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr2.csv")) %>%
  select(all_of(keep))
d <- read.csv(paste0("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr3.csv")) %>%
  select(all_of(keep))
e <- read.csv(paste0("data/invitro-vs-mc_results_D7_down-sampled_expanded_gr4.csv")) %>%
  select(all_of(keep))

a <- a %>% mutate(Clone.Size = ">0")
b <- b %>% mutate(Clone.Size = ">1")
c <- c %>% mutate(Clone.Size = ">2")
d <- d %>% mutate(Clone.Size = ">3")
e <- e %>% mutate(Clone.Size = ">4")

summary <- rbind(a,b,c,d,e)

beta <- summary %>% filter(chain == "beta") %>%
  pivot_longer(cols = starts_with("pct"),
               names_to = "TCR.class",
               values_to = "percentage") %>%
  mutate(TCR.class = recode(TCR.class,
                            pct.PPD.only = "private",
                            pct.mc = "public"))

# plot (save as svg 600x350)
ggplot(beta, aes(x=TCR.class,y=percentage,colour=TCR.class)) +
  geom_line(aes(group=UIN),color="black")+
  geom_point(size=2,position = position_dodge(width=0.7)) +
  facet_wrap(~Clone.Size,scales="free_x",strip.position = "bottom",ncol=5)+
  My_Theme +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0,units="line"))+
  scale_colour_manual(values = c("orange","navy")) +
  labs(y="% of total CDR3s (down-sampled)",
       x = "Clone size") +
  ggtitle("Day 7 TST")
