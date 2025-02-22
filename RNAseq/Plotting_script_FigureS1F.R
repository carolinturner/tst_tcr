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
dat <- read.csv("pathway_analysis.csv")
# retain for each pathway group (= cluster cut) the name and statistics of the biggest annotated pathway
dat <- dat %>% filter(!is.na(nAnno))

result <- dat %>% 
  group_by(group) %>%
  filter(nAnno == max(nAnno))

result$Group = factor(result$Group, levels = c("Day 2 TST","Day 7 TST")) # change levels according to input data
rgb.palette1 <- colorRampPalette(c("blue", "red"), space = "rgb")

ggplot() +
  geom_point(data=result,aes(x=Group,y=reorder(name,nAnno),size=log2(nOverlap),colour=zscore)) +
  My_Theme+
  theme(axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = -45, hjust = 0)) +
  xlab("") + ylab("") +
  scale_colour_gradientn(colours=rgb.palette1(20)) +
  scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE)
