library(tidyverse)
library(epitools)
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
  panel.spacing = unit(0.05,"cm")
  )

# load data and calculate odds ratios
d1 <- read.csv("data/Poisson_expanded_D7-v-D2_alpha.csv")
ppd <- d1 %>% filter(PPD == PPD) %>% nrow()
tt <- d1 %>% filter(TT == TT) %>% nrow()
M1 <- matrix(c(ppd,nrow(d1)-ppd,tt,nrow(d1)-tt),ncol=2)
OR1 <- as.data.frame(t(round(oddsratio(M1)$measure[2,],digits = 1)))
dat1 <- data.frame(dataset = rep("expanded_alpha",2),
                   Stimulant = c("PPD","TT"),
                   Value = c(ppd/nrow(d1)*100,
                             tt/nrow(d1)*100))

d2 <- read.csv("data/Poisson_non-expanded_D7-v-D2_alpha.csv")
ppd <- d2 %>% filter(PPD == PPD) %>% nrow()
tt <- d2 %>% filter(TT == TT) %>% nrow()
M2 <- matrix(c(ppd,nrow(d2)-ppd,tt,nrow(d2)-tt),ncol=2)
OR2 <- as.data.frame(t(round(oddsratio(M2)$measure[2,],digits = 1)))
dat2 <- data.frame(dataset = rep("non-expanded_alpha",2),
                   Stimulant = c("PPD","TT"),
                   Value = c(ppd/nrow(d2)*100,
                             tt/nrow(d2)*100))

d3 <- read.csv("data/Poisson_expanded_D7-v-D2_beta.csv")
ppd <- d3 %>% filter(PPD == PPD) %>% nrow()
tt <- d3 %>% filter(TT == TT) %>% nrow()
M3 <- matrix(c(ppd,nrow(d3)-ppd,tt,nrow(d3)-tt),ncol=2)
OR3 <- as.data.frame(t(round(oddsratio(M3)$measure[2,],digits = 1)))
dat3 <- data.frame(dataset = rep("expanded_beta",2),
                   Stimulant = c("PPD","TT"),
                   Value = c(ppd/nrow(d3)*100,
                             tt/nrow(d3)*100))

d4 <- read.csv("data/Poisson_non-expanded_D7-v-D2_beta.csv")
ppd <- d4 %>% filter(PPD == PPD) %>% nrow()
tt <- d4 %>% filter(TT == TT) %>% nrow()
M4 <- matrix(c(ppd,nrow(d4)-ppd,tt,nrow(d4)-tt),ncol=2)
OR4 <- as.data.frame(t(round(oddsratio(M4)$measure[2,],digits = 1)))
dat4 <- data.frame(dataset = rep("non-expanded_beta",2),
                   Stimulant = c("PPD","TT"),
                   Value = c(ppd/nrow(d4)*100,
                             tt/nrow(d4)*100))

# summary
d <- rbind(dat1,dat2,dat3,dat4)
OR <- rbind(OR1,OR2,OR3,OR4)
rownames(OR) <- c("expanded_alpha","non-expanded_alpha","expanded_beta","non-expanded_beta")
write.csv(OR,"data/FigureS8B_odds-ratio.csv")

# prep for plotting
df <- d %>%
  mutate(chain = c(rep("alpha",4),rep("beta",4)),
         dataset = recode(dataset,
                          'expanded_alpha' = "expanded CDR3s",
                          'non-expanded_alpha' = "non-expanded CDR3s",
                          'expanded_beta' = "expanded CDR3s",
                          'non-expanded_beta' = "non-expanded CDR3s")) 
  
dat_text <- data.frame(
  Value = c(rep(70,4)),
  Stimulant = c(rep("TT",4)),
  chain = c(rep("alpha",2),rep("beta",2)),
  dataset = c(rep(c("expanded CDR3s","non-expanded CDR3s"),2)),
  label = c(paste0("OR:\n",OR[1,1]," (",OR[1,2],"-",OR[1,3],")"),
            paste0("OR:\n",OR[2,1]," (",OR[2,2],"-",OR[2,3],")"),
            paste0("OR:\n",OR[3,1]," (",OR[3,2],"-",OR[3,3],")"),
            paste0("OR:\n",OR[4,1]," (",OR[4,2],"-",OR[4,3],")"))
)

# plot
pS8B <- ggplot(df, aes(x=Stimulant,y=Value))+
  geom_col()+
  facet_grid(chain ~ dataset)+
  geom_text(data=dat_text,label = dat_text$label,size = 6/.pt)+
  My_Theme+
  labs(y="Percentage of antigen-reactive CDR3s")

ggarrange(pS8B,
          labels = c("B"),
          font.label = list(size = 10, face = "bold", colour = "black"))
ggsave("figures/FigureS8B.svg", 
       units = "cm", width = 8, height =8 , dpi=300)
