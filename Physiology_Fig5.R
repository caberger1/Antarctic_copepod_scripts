#Make multipanel plot of C. acutus and C. propinquus physiology
# 6 panels: FIG (food in gut), reprod, CS (citrate synthase)

rm(list = ls())
library("ggplot2")
library("dplyr")
library(ggpubr)
library(export)
setwd("yourpath")
#throughout script replace "yourpath" with path to file


#Food in gut ("FIG")
data <- read.csv("Acutus_Site200_FIG_Table.csv")
data <- data[,2:4]

data$LatPhoto = factor(data$LatPhoto) # convert to nominal factor
data$SiteGroup = factor(data$SiteGroup, ordered = TRUE,levels = c("S200", "Other")) 
levels(data$SiteGroup) <- c("200.xxx", "Others")
data$Final_FIG = factor(data$Final_FIG, ordered = TRUE, levels = c("Yes", "No")) 
summary(data)

df <- data %>%
  group_by(SiteGroup, Final_FIG, .drop= FALSE) %>%
  summarise(counts = n())  

legend_title <- "Food in Gut"
p <- ggplot(df, aes(x = SiteGroup, y = counts)) +
  geom_bar(
    aes(color = Final_FIG, fill = Final_FIG),
    stat = "identity", position = position_dodge(0.8),
    width = 0.7
  ) +
  xlab("Station Group")+
  ylab("# Individuals")+
  theme_classic()+ 
  scale_fill_manual(legend_title, values = c("#228833",  "#BBBBBB"))+
  scale_color_manual(values = c("black",  "black"))
p <- p + guides(color = "none")

# Reproduction
library(RColorBrewer)
Repro_data <- read.csv(file= "C_acutus_reprod_summarized.csv")
Repro_data$Stage = factor(Repro_data$Stage, ordered = TRUE,levels = c("1", "2", "3")) 
Repro_data$Site = factor(Repro_data$Site, ordered = TRUE,levels = c("Site200", "Others")) 
levels(Repro_data$Site) <- c("200.xxx", "Others")
legend_title <- "Reproductive Stage"

q <- ggplot(Repro_data, aes(x = Site, y = Number)) +
  geom_bar(
    aes(color = Stage, fill = Stage),
    stat = "identity")+
  xlab("Station Group")+
  ylab("# Individuals")+
  scale_fill_brewer(legend_title, palette="Reds")+
  scale_color_manual(values = c("black",  "black", "black"))+
  theme_classic()
q <- q + guides(color = "none")

  
  # Citrate synthase
  Ca_Station <- read.csv("CS_Ca_Station.csv", header=TRUE)
  Ca_Station <- Ca_Station[which (Ca_Station$Cryo != "LC9"),] #remove sample with a Gigas in it
  Ca_Station$Station = factor(Ca_Station$Station)
  Ca_Station$Group = factor(Ca_Station$Group)
  levels(Ca_Station$Station) <- c("600.100", "200.040", "200.-040", "200.000", "000.100", "-100.100", "500.100", "117")
  CaSeq <- c("200.040", "200.000", "000.100", "-100.100") # 4 RNA-seq sites with at least 3 (4) CS measurements
  vec <- c('A','a','C')
  
  Ca_SeqSta <- Ca_Station[Ca_Station$Station %in% CaSeq,]
  
 r<- ggplot(Ca_SeqSta, aes(x=Station, y=CS_weight, fill=Station)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#CCBB44", "#228833", "#4477AA", "#BBBBBB"))+
    geom_point(color="black", size=0.9, alpha=0.9) +
    xlab("Station")+
    ylab("CS, U/mg")+
    theme_classic()+
    theme(legend.position="none") 
 
 # make a version collapsed by station group
 r2<- ggplot(Ca_SeqSta, aes(x=Group, y=CS_weight, fill=Group)) +
   scale_fill_manual(values = c("#228833", "#BBBBBB")) +
   geom_boxplot() +
   geom_point(color="black", size=0.9, alpha=0.9) +
   xlab("Station Group")+
   ylab("CS, U/mg")+
   theme_classic()+
   theme(legend.position="none") 
 
#Propinquus
 data <- read.csv("CpFIG.csv")
 #data <- data[,2:4]
 
 data$LatPhoto = factor(data$LatPhoto) # convert to nominal factor
 data$Station = factor(data$Station, ordered=TRUE, levels=c("200", "100.18", "0.1", "-100.1")) 
 #data$SiteGroup = factor(data$SiteGroup, ordered = TRUE,levels = c("S200", "Other")) 
 #levels(data$SiteGroup) <- c("200.000 & 200.040", "Other Sites")
 data$Final_FIG = factor(data$Final_FIG, ordered = TRUE, levels = c("Yes", "No")) 
 summary(data)
 
 df <- data %>%
   group_by(Station, Final_FIG, .drop= FALSE) %>%
   summarise(counts = n())  
 
 legend_title <- "Food in Gut"
 s <- ggplot(df, aes(x = Station, y = counts)) +
   geom_bar(
     aes(color = Final_FIG, fill = Final_FIG),
     stat = "identity", position = position_dodge(0.8),
     width = 0.7
   ) +
   xlab("Station")+
   ylab("# Individuals")+
   theme_classic()+ 
   theme(legend.position="none")+ 
   scale_fill_manual(legend_title, values = c("#228833",  "#BBBBBB"))+
   scale_color_manual(values = c("black",  "black")) 
 s <- s + guides(color = "none")
 
 # Reproduction
 library(RColorBrewer)
 Repro_data <- read.csv("C_propinquus_reprod_summarized.csv")
 
 
 Repro_data$Stage = factor(Repro_data$Stage, ordered = TRUE,levels = c("0", "1", "2", "3")) 
 Repro_data$Station = factor(Repro_data$Site, ordered = TRUE,levels = c("200", "100.18", "0.1", "-100.1")) 
 #levels(Repro_data$Site) <- c("200.000 & 200.040", "Other Sites")
 legend_title <- "Reproductive Stage"
 
 t <- ggplot(Repro_data, aes(x = Station, y = Number)) +
   geom_bar(
     aes(color = Stage, fill = Stage),
     stat = "identity")+
   xlab("Station")+
   ylab("# Individuals")+
   scale_fill_brewer(legend_title, palette="Reds")+
   scale_color_manual(values = c("black",  "black", "black"))+
   theme_classic()+
   theme(legend.position="none") 
 t <- t + guides(color = "none") 
 

 Cp_Station <- read.csv("CS_Cp_Station.csv", header=TRUE)
 Cp_Station$Station = factor(Cp_Station$Station)
 levels(Cp_Station$Station) <- c("100.180", "000.100", "-100.100", "-100.160")
  Cp_3Stations <- Cp_Station[ which(Cp_Station$Station!='100.180'),]
 Cp_2Stations<- Cp_3Stations[ which(Cp_3Stations$Station!='-100.160'),] 
 
 u <- Cp_2Stations %>%
   ggplot( aes(x=Station, y=CS_weight, fill=Station)) +
   geom_boxplot() +
   scale_fill_manual(values = c("#4477AA", "#BBBBBB")) +
   geom_point(color="black", size=0.9, alpha=0.9) +
   xlab("Station")+
   ylab("CS, U/mg")+
   ylim(1,5)+
   theme_classic()+
   theme(legend.position="none")
   
 
 #Plot as desired
 #in powerpoint
 ggarrange(p, s, q, t, r2, u,
           labels = c("A", "B", "C", "D", "E", "F"),
           ncol = 2, nrow = 3)
 #graph2ppt(file="test.pptx", width=8, height=7) 
 
 #as png
 #png(filename = "C:/Users/atarrant/Documents/Paper_Ant_Spatial/Figures_and_Tables/Compbined_physiology.png")
 ggarrange(p, s, q, t, r2, u,
            labels = c("A", "B", "C", "D", "E", "F"),
            ncol = 2, nrow = 3)
  #dev.off()
 
 #as pdf
  #pdf(file = "C:/Users/atarrant/Documents/Paper_Ant_Spatial/Figures_and_Tables/Combined_physiology.pdf")
  ggarrange(p, s, q, t, r2, u,
            labels = c("A", "B", "C", "D", "E", "F"),
            ncol = 2, nrow = 3)
  #dev.off()