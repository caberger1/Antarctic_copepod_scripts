#Script makes the heatmap (Fig 3 in original ms) showing correlations among environmental parameters and copepod abundances
#throughout script replace "yourpath" with path to file

library(dplyr)
library("Hmisc")
library(corrplot)
library(ggplot2) 

rm(list = ls())
setwd("yourpath")

data <- read.csv("CorrPlot_factors.csv", header = TRUE, stringsAsFactors = TRUE) %>% 
 select(Calanoides, Calanus, Copepod_vol, Temp_avg, Chlorophyll_int, Latitude, Depth) %>% 
  rename(Chlorophyll=Chlorophyll_int, Temperature = Temp_avg, Copepods = Copepod_vol)

data2 <- data %>% 
  select(Calanoides, Calanus, Latitude, Chlorophyll, Temperature)

#Correlations among environmental variables

res3 <- rcorr(as.matrix(data),type = "spearman")
diag(res3$r) = NA
col1<- colorRampPalette(c("blue3", "white", "red3"))

#output as desired
#pdf(file = "Figures_and_Tables/CorrPlot.pdf")
corrplot(res3$r, type = "upper", col = col1(100),
             tl.col = "black", tl.srt = 45, na.label="-", p.mat = res3$P, sig.level = .05)
#dev.off()

res4 <- rcorr(as.matrix(data2),type = "spearman")
diag(res4$r) = NA
col1<- colorRampPalette(c("blue3", "white", "red3"))
z<-corrplot(res4$r, type = "upper", col = col1(100),
            tl.col = "black", tl.srt = 45, na.label="-", p.mat = res4$P, sig.level = .05)


