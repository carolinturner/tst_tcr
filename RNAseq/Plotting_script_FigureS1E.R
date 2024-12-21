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

# read in data
dat <- read.table("SARtools_output_D7vsD2/D7vsD2.complete.txt", header = T, row.names = 1)

# add a column to specify if genes are up- or down-regulated (using as thresholds: log2 FC 1 and padj 0.05)
dat$Genes <- "no differential expression"
dat$Genes[dat$log2FoldChange >= 1 & dat$padj < 0.05] <- "enriched in Day 7 TST"
dat$Genes[dat$log2FoldChange <= -1 & dat$padj < 0.05] <- "enriched in Day 2 TST"

# check number of up- and down-regulated genes
sum(dat$Genes == "enriched in Day 7 TST")
sum(dat$Genes == "enriched in Day 2 TST")

# make a custom colour vector
mycolour <- c("blue","red","black")
names(mycolour) <- c("enriched in Day 2 TST","enriched in Day 7 TST","no differential expression")

# make "volcano plot" (= scatter plot with log2 fold change on x-axis and adj p value on y-axis)
ggplot(dat, aes(x=log2FoldChange,y=-log10(padj),col=Genes)) +
  geom_point(alpha=0.5) +
  geom_vline(xintercept = c(-1,1), col="grey") +
  geom_hline(yintercept = -log10(0.05), col="grey") +
  My_Theme +
  scale_color_manual(values = mycolour)  
