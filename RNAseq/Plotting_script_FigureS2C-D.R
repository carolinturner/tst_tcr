library(tidyverse)
library(reshape2)

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

# read in module matrix and metadata
dat <- read.csv("data/TST_fdrsig_UR_module_expression_per_sample.csv", header=T, row.names = 1)
meta <- read.csv("data/RNAseq_metadata.csv", header=T, row.names = 1) %>%
  select(Stimulant)

### Mann Whitney tests with FDR correction
meta.d2 <- subset(meta, Stimulant == "TST_D2")
meta.d7 <- subset(meta, Stimulant == "TST_D7")

dat.d2 <- dat[,row.names(meta.d2)]
dat.d7 <- dat[,row.names(meta.d7)]

# Statistical comparison D2 vs D7
wc <- sapply(1:nrow(dat), function(i){
  wilcox.test(as.numeric(dat.d2[i,]),
              as.numeric(dat.d7[i,]))$p.value})  

summary <- data.frame(row.names(dat),wc)
colnames(summary) <- c("module","MW.p")

# multiple comparison correction
summary$BH.fdr <- p.adjust(summary$MW.p, method = "BH")

# Calculate fold change between median values for each module and add to summary
summary$median.d2 <- apply(dat.d2, 1, median)
summary$median.d7 <- apply(dat.d7, 1, median)
summary$LFC.d2vsd7 <- summary$median.d2 - summary$median.d7

# set module name as rowname
summary <- summary %>% column_to_rownames("module") 

# add grouping based on significance and direction of change
summary$group <- ifelse(summary$BH.fdr > 0.05, "non.sig",
                        ifelse(summary$LFC.d2vsd7 < 0, "up.d7", "up.d2"))

# reformat for plotting
d7 <- summary %>%
  filter(BH.fdr < 0.05 & LFC.d2vsd7 < 0) %>%
  select(median.d2, median.d7) %>%
  rename("Day 2 TST" = median.d2, "Day 7 TST" = median.d7) %>%
  rownames_to_column("module") %>%
  pivot_longer(cols = c("Day 2 TST","Day 7 TST"),names_to = "Sample")
d2 <- summary %>%
  filter(BH.fdr < 0.05 & LFC.d2vsd7 > 0) %>%
  select(median.d2, median.d7) %>%
  rename("Day 2 TST" = median.d2, "Day 7 TST" = median.d7) %>%
  rownames_to_column("module") %>%
  pivot_longer(cols = c("Day 2 TST","Day 7 TST"),names_to = "Sample")

# plots
ggplot(d2, aes(x=Sample, y=value, group=module),colour="black") +
  geom_line(linewidth =0.7) +
  geom_point(size =2,alpha=0.5) +
  ylab("Median module TPM") +
  My_Theme+
  theme(legend.position="none")

ggplot(d7, aes(x=Sample, y=value, group=module, color=module)) +
  geom_line(linewidth =1) +
  geom_point(size =2) +
  ylab("Median module TPM") +
  My_Theme +
  labs(colour = "Upstream regulator module") 
