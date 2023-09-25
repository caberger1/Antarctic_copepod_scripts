#Script plots depth profiles of environmental data for selected stations, as shown in Figure 2 of the original manuscript
#throughout script replace "yourpath" with path to file

library(ggplot2)
library(ggpubr)

# Read in data. 
rm(list = ls())
setwd("yourpath")
#"Sta numbers" are CTD cast numbers. 
#Sta16=400.040 Sta22=200.040 Sta28=200.000 Sta34=100.180 Sta36=000.100 Sta38=-100.000 Sta48=616.040
Sta16 <- read.table("Sta16_edit.txt", header = TRUE)
Sta22 <- read.table("Sta22_edit.txt", header = TRUE)
Sta28 <- read.table("Sta28_edit.txt", header = TRUE)
Sta34 <- read.table("Sta34_edit.txt", header = TRUE)
Sta36 <- read.table("Sta36_edit.txt", header = TRUE)
Sta38 <- read.table("Sta38_edit.txt", header = TRUE)
Sta48 <- read.table("Sta48_edit.txt", header = TRUE)

a<-ggplot()+
  scale_y_reverse(limits=c(200,0))+
  ylab("Depth, m")+
  xlab("Temperature, \u00B0C")+
  labs(colour="Site")+
  geom_path(data=Sta48, aes(x=T090C, y=DepSM), color = "#AA3377", size=1.5)+
  geom_path(data=Sta16, aes(x=T090C, y=DepSM), color = "#66CCEE", size=1.5)+
  geom_path(data=Sta22, aes(x=T090C, y=DepSM), color = "#CCBB44",  size=1.5)+
  geom_path(data=Sta28, aes(x=T090C, y=DepSM), color = "#228833", size=1.5)+
  geom_path(data=Sta34, aes(x=T090C, y=DepSM), color = "#EE6677", size=1.5)+
  geom_path(data=Sta36, aes(x=T090C, y=DepSM), color = "#4477AA", size=1.5)+
  geom_path(data=Sta38, aes(x=T090C, y=DepSM), color = "#707070", size=1.5)+
  theme_classic()+
  theme(axis.text.x=element_text(size=rel(1.5)), axis.title.x=element_text(size=rel(1.5)))+
  theme(axis.text.y=element_text(size=rel(1.5)), axis.title.y=element_text(size=rel(1.5)))+
  theme(legend.position="none")
  theme(legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.5)))+
  
b<-ggplot()+
  scale_y_reverse(limits=c(200,0))+
  ylab("Depth, m")+
  xlab(expression("Chlorophyll mg/"~m^3))+
  labs(colour="Site")+
  geom_path(data=Sta48, aes(x=FlECO.AFL, y=DepSM, color=" 616.040"), size=1.5)+
  geom_path(data=Sta16, aes(x=FlECO.AFL, y=DepSM, color=" 400.040"), size=1.5)+
  geom_path(data=Sta22, aes(x=FlECO.AFL, y=DepSM, color=" 200.040"),  size=1.5)+
  geom_path(data=Sta28, aes(x=FlECO.AFL, y=DepSM, color=" 200.000"), size=1.5)+
  geom_path(data=Sta34, aes(x=FlECO.AFL, y=DepSM, color=" 100.180"), size=1.5)+
  geom_path(data=Sta36, aes(x=FlECO.AFL, y=DepSM, color=" 000.100"), size=1.5)+
  geom_path(data=Sta38, aes(x=FlECO.AFL, y=DepSM, color="-100.000"), size=1.5)+
  theme_classic()+
  theme(axis.text.x=element_text(size=rel(1.5)), axis.title.x=element_text(size=rel(1.5)))+
  theme(axis.text.y=element_text(size=rel(1.5)), axis.title.y=element_text(size=rel(1.5)))+
  theme(legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.5)))+
  theme(legend.position = c(0.8,0.4))+
  scale_color_manual(values = c("#AA3377", "#66CCEE", "#CCBB44", "#228833", "#EE6677", "#4477AA", "#BBBBBB"),
                     breaks = c(" 616.040", " 400.040", " 200.040", " 200.000", " 100.180", " 000.100", "-100.100"))

c<-ggplot()+
  scale_y_reverse(limits=c(200,0))+
  ylab("Depth, m")+
  xlab("Salinity")+
  labs(colour="Site")+
  geom_path(data=Sta48, aes(x=Sal00.1, y=DepSM, color=" 616.040"), size=1.5)+
  geom_path(data=Sta16, aes(x=Sal00.1, y=DepSM, color=" 400.040"), size=1.5)+
  geom_path(data=Sta22, aes(x=Sal00.1, y=DepSM, color=" 200.040"),  size=1.5)+
  geom_path(data=Sta28, aes(x=Sal00.1, y=DepSM, color=" 200.000"), size=1.5)+
  geom_path(data=Sta34, aes(x=Sal00.1, y=DepSM, color=" 100.180"), size=1.5)+
  geom_path(data=Sta36, aes(x=Sal00.1, y=DepSM, color=" 000.100"), size=1.5)+
  geom_path(data=Sta38, aes(x=Sal00.1, y=DepSM, color="-100.000"), size=1.5)+
  theme_classic()+
  theme(axis.text.x=element_text(size=rel(1.5)), axis.title.x=element_text(size=rel(1.5)))+
  theme(axis.text.y=element_text(size=rel(1.5)), axis.title.y=element_text(size=rel(1.5)))+
  theme(legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.5)))+
  theme(legend.position = c("none"))+
  scale_color_manual(values = c("#AA3377", "#66CCEE", "#CCBB44", "#228833", "#EE6677", "#4477AA", "#BBBBBB"),
                     breaks = c(" 616.040", " 400.040", " 200.040", " 200.000", " 100.180", " 000.100", "-100.100"))

ggarrange(a,b,c, ncol=2, nrow=2, labels = c("a", "b", "c"), font.label = list(size=20, face = "plain"))

#write out file as desired
#pdf(file="Profiles.pdf")
ggarrange(a,b, ncol=2, nrow=1, labels = c("A", "B"), font.label = list(size=20, face = "plain"))
#dev.off()
