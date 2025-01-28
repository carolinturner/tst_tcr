library(tidyverse)
library(epitools)

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

# load data and calculate odds ratios
d1 <- read.csv("data/Poisson_expanded_D7-v-D2_alpha.csv")
ppd <- d1 %>% filter(PPD == PPD) %>% nrow()
tt <- d1 %>% filter(TT == TT) %>% nrow()
M1 <- matrix(c(ppd,nrow(d1)-ppd,tt,nrow(d1)-tt),ncol=2)
OR1 <- as.data.frame(t(round(oddsratio(M1)$measure[2,],digits = 2)))
dat1 <- data.frame(dataset = rep("expanded_alpha",2),
                   Stimulant = c("PPD","TT"),
                   Value = c(ppd/nrow(d1)*100,
                             tt/nrow(d1)*100))

d2 <- read.csv("data/Poisson_non-expanded_D7-v-D2_alpha.csv")
ppd <- d2 %>% filter(PPD == PPD) %>% nrow()
tt <- d2 %>% filter(TT == TT) %>% nrow()
M2 <- matrix(c(ppd,nrow(d2)-ppd,tt,nrow(d2)-tt),ncol=2)
OR2 <- as.data.frame(t(round(oddsratio(M2)$measure[2,],digits = 2)))
dat2 <- data.frame(dataset = rep("non-expanded_alpha",2),
                   Stimulant = c("PPD","TT"),
                   Value = c(ppd/nrow(d2)*100,
                             tt/nrow(d2)*100))

d3 <- read.csv("data/Poisson_expanded_D7-v-D2_beta.csv")
ppd <- d3 %>% filter(PPD == PPD) %>% nrow()
tt <- d3 %>% filter(TT == TT) %>% nrow()
M3 <- matrix(c(ppd,nrow(d3)-ppd,tt,nrow(d3)-tt),ncol=2)
OR3 <- as.data.frame(t(round(oddsratio(M3)$measure[2,],digits = 2)))
dat3 <- data.frame(dataset = rep("expanded_beta",2),
                   Stimulant = c("PPD","TT"),
                   Value = c(ppd/nrow(d3)*100,
                             tt/nrow(d3)*100))

d4 <- read.csv("data/Poisson_non-expanded_D7-v-D2_beta.csv")
ppd <- d4 %>% filter(PPD == PPD) %>% nrow()
tt <- d4 %>% filter(TT == TT) %>% nrow()
M4 <- matrix(c(ppd,nrow(d4)-ppd,tt,nrow(d4)-tt),ncol=2)
OR4 <- as.data.frame(t(round(oddsratio(M4)$measure[2,],digits = 2)))
dat4 <- data.frame(dataset = rep("non-expanded_beta",2),
                   Stimulant = c("PPD","TT"),
                   Value = c(ppd/nrow(d4)*100,
                             tt/nrow(d4)*100))

# summary
d <- rbind(dat1,dat2,dat3,dat4)
OR <- rbind(OR1,OR2,OR3,OR4)
rownames(OR) <- c("expanded_alpha","non-expanded_alpha","expanded_beta","non-expanded_beta")
write.csv(OR,"data/FigureS6_odds-ratio.csv")

# plot (save as svg 400x400)
df <- d %>%
  mutate(chain = c("alpha","alpha","alpha","alpha","beta","beta","beta","beta"),
         dataset = recode(dataset,
                          'expanded_alpha' = "expanded TCRs",
                          'non-expanded_alpha' = "non-expanded TCRs",
                          'expanded_beta' = "expanded TCRs",
                          'non-expanded_beta' = "non-expanded TCRs"))

ggplot(df, aes(x=Stimulant,y=Value))+
  geom_col()+
  facet_grid(chain ~ dataset)+
  My_Theme+
  labs(y="Percentage of antigen-reactive TCRs")
