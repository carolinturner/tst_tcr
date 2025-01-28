# TCRseq_script_7: Expansion of Mtb-reactive TCRs in TST between day 2 and day 7
# Script adapted from Benny Chain

library(data.table)
library(tidyverse)
library(mgcv)
#library(textshape)

# load meta data and extract UINs with paired TST samples
meta <- read.csv("data/TCRseq_metadata.csv")
d2 <- meta %>% filter(tissue == "TST_D2") %>% pull(UIN) %>% unique()
d7 <- meta %>% filter(tissue == "TST_D7") %>% pull(UIN) %>% unique()
paired <- intersect(d2,d7)

spl <- meta %>%
  filter(tissue %in% c("TST_D2","TST_D7") & UIN %in% paired) %>%
  pull(UIN) %>% unique()

# load data and subset to samples of relevance
dat.a <- fread("data/combined_alpha.csv.gz") %>%
  filter(tissue %in% c("TST_D2","TST_D7") & UIN %in% spl)
dat.b <- fread("data/combined_beta.csv.gz") %>%
  filter(tissue %in% c("TST_D2","TST_D7") & UIN %in% spl)

# loop through each chain
for (chain in c("alpha","beta")){
  dat <- if(chain == "alpha"){dat.a}else{dat.b}
  
  # make summary dataframes
  exp <- data.frame(CDR3=character(),subject=character())
  non.exp <- data.frame(CDR3=character(),subject=character())
  
  # load PPD- and TT-reactive CDR3s
  PPD.CDR3 <- as.character(read.csv(paste0("data/PPD_",chain,".csv"), header = F)$V1)
  TT.CDR3 <- as.character(read.csv(paste0("data/TT_",chain,".csv"), header = F)$V1)
  
  # loop through each subject
  for(i in spl){
    # select data
    D2 <- dat %>% filter(UIN == i & tissue == "TST_D2")
    D7 <- dat %>% filter(UIN == i & tissue == "TST_D7")
    # sum up number of CDR3s independent of gene rearrangement
    D2 <- D2 %>% group_by(junction_aa) %>% summarise(total = sum(duplicate_count)) %>% ungroup()
    D7 <- D7 %>% group_by(junction_aa) %>% summarise(total = sum(duplicate_count)) %>% ungroup()
    # rename columns
    colnames(D2) <- c("CDR3", "count_D2")
    colnames(D7) <- c("CDR3", "count_D7")
    # merge data, keeping all clones from both samples
    merge <- merge(D7,D2,all=T)
    
    ### Poisson distribution to identify expansion between D2 and D7 - adapted from Benny Chain
    
    # replace missing clones with median count of the respective samples (usually =1)
    merge$count_D2[is.na(merge$count_D2)] <- median(D2$count_D2)
    merge$count_D7[is.na(merge$count_D7)] <- median(D7$count_D7)
    # add PPD reactivity
    merge$PPD <- ifelse(merge$CDR3 %in% PPD.CDR3, "PPD", "NA")
    
    # normalise data (counts/million)
    cpm <- merge
    cpm$count_D2 <- cpm$count_D2/sum(cpm$count_D2)*1000000
    cpm$count_D7 <- cpm$count_D7/sum(cpm$count_D7)*1000000
    # log2 transform data
    log2 <- cpm
    log2$count_D2 <- log2(log2$count_D2)
    log2$count_D7 <- log2(log2$count_D7)
    
    # calculate Poisson distribution
    #probability cut-off for the Poisson set to 0.0001 because it gave very few TCRs in control
    prob<-0.0001
    #calculate the TCR count that would give a significant value compared to 1,2,4,8,16,32,64.
    poiss_l<- qpois(prob,c(1,2^(0:8)),lower.tail=FALSE) 
    
    #error margins - calculating the confidence intervals that, if you observe a TCR m times on D2,
    #you will observe it n times on D7, if the only process at work is random Poisson sampling.
    #multiply by a scalar because all the raw counts are normalised to TCR counts/million.
    
    # scaling constants - I have been using the maximum scaling factor for the two samples being compared.
    # but not sure whether we could use different scaling in different direction
    
    scalar_2<-1E6/length(D2$CDR3)
    scalar_1<-1E6/length(D7$CDR3)
    scalar<-max(scalar_1,scalar_2)
    
    # plot data
    png(paste0("data/Poisson_plot_",i,"_",chain,".png"),units="in", width=5, height=5,res=300)
    plot(log2$count_D2,log2$count_D7,xlim = c(2,16),ylim=c(2,16),
         ylab="log2 TCRs/million (TST_D7)",xlab="log2 TCRs/million (TST_D2)",
         col=factor(log2$PPD))
    # add the error limits as calculated above using Poisson distribution
    #lower limit
    points(log2(c(min(log2$count_D2),(2^c(0:8,8))*1E6/sum(D2$count_D2))),log2(c(poiss_l*scalar,2^16)),col="blue",pch=19,type="l", lty = 2,lwd=2)
    
    #upper limits
    points(log2(c(poiss_l*scalar,2^16)),log2(c(min(log2$count_D7),(2^c(0:8,8))*1E6/sum(D7$count_D7))),col="blue",pch=19,type="l", lty = 2,lwd=2)
    
    dev.off()
    
    # extract numbers
    bnd_up<-cbind(log2(c(min(log2$count_D2),(2^c(0:8,8))*1E6/sum(D2$count_D2),min(log2$count_D2))),log2(c(poiss_l*scalar,2^16,2^16)))
    
    i_up<-which(in.out(bnd_up,cbind(log2$count_D2,log2$count_D7)))
    
    #collect all expanded TCRs
    TCR_up<-as.data.frame(cbind(log2$CDR3[i_up]))
    TCR_up$UIN <- i
    colnames(TCR_up)[1] <- "CDR3"
    
    exp <- rbind(exp,TCR_up)
    
    # collect all non-expanded TCRs
    TCR_ne <- as.data.frame(cbind(log2$CDR3[-i_up]))
    TCR_ne$UIN <- i
    colnames(TCR_ne)[1] <- "CDR3"
    non.exp <- rbind(non.exp,TCR_ne)
  }
  
  # convert expanded TCRs to matrix and annotate with Ag reactivity
  wide <- reshape(exp, idvar="CDR3",v.names = "UIN", timevar = "UIN", direction = "wide")
  
  wide <- wide %>% remove_rownames() %>% column_to_rownames("CDR3")
  wide[!is.na(wide)] <-1
  wide[is.na(wide)] <- 0
  
  wide <- wide %>% tibble::rownames_to_column("CDR3")
  wide$PPD <- ifelse(wide$CDR3 %in% PPD.CDR3, "PPD", "NA")
  wide$TT <- ifelse(wide$CDR3 %in% TT.CDR3, "TT", "NA")
  
  write.csv(wide,paste0("data/Poisson_expanded_D7-v-D2_",chain,".csv"),row.names = F)
  
  # convert non-expanded TCRs to matrix and annotate with Ag reactivity
  wide.ne <- reshape(non.exp, idvar="CDR3",v.names = "UIN", timevar = "UIN", direction = "wide")
  
  wide.ne <- wide.ne %>% remove_rownames() %>% column_to_rownames("CDR3")
  wide.ne[!is.na(wide.ne)] <-1
  wide.ne[is.na(wide.ne)] <- 0
  
  wide.ne <- wide.ne %>% tibble::rownames_to_column("CDR3")
  wide.ne$PPD <- ifelse(wide.ne$CDR3 %in% PPD.CDR3, "PPD", "NA")
  wide.ne$TT <- ifelse(wide.ne$CDR3 %in% TT.CDR3, "TT", "NA")
  
  write.csv(wide.ne,paste0("data/Poisson_non-expanded_D7-v-D2_",chain,".csv"),row.names = F)
}
