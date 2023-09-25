#Script makes the subplots (Principal components analysis) shown in Fig. 4. Subplots were assembled into figure outside of R. 
library(ggplot2)

#replace "yourpath" with path to file
setwd("yourpath")
df <- read.csv(file="Acutus_PC_data.csv")
df$Station <- factor(df$Station, ordered = TRUE,levels = c("616.04", "400.04", "200.04", "200", "0.1", "-100.1"))
pdf(file="Acutus_PCA.pdf")
ggplot(df) + 
  geom_point(aes(x=PC1, y=PC2, color = Station), size=5)+ 
  xlab("PC1 (13.3%)")+
  ylab("PC2 (7.3%)")+
  theme_classic()+ 
  theme(legend.position = "none")+
  theme(aspect.ratio=3/6)+
  theme(axis.text.x=element_text(size=rel(1.5)), axis.title.x=element_text(size=rel(1.5)))+
  theme(axis.text.y=element_text(size=rel(1.5)), axis.title.y=element_text(size=rel(1.5)))+
  theme(legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.2)))+
  scale_color_manual(values = c("#AA3377", "#66CCEE", "#CCBB44", "#228833", "#4477AA", "#BBBBBB"))
dev.off()

df <- read.csv(file = "Cp_PCA_sites.csv")
df$Station <- factor(df$Station, ordered = TRUE,levels = c("-100.1", "0.1", "100.18", "200"))
pdf(file="Propinquus_PCA.pdf")
ggplot(df) + 
  geom_point(aes(x=PC1, y=PC2, color = Station), size=5)+ 
  xlab("PC1 (9.7%)")+
  ylab("PC2 (9.0%")+
  theme_classic()+ 
  theme(legend.position = "none")+
  theme(aspect.ratio = 3/6)+
  theme(axis.text.x=element_text(size=rel(1.5)), axis.title.x=element_text(size=rel(1.5)))+
  theme(axis.text.y=element_text(size=rel(1.5)), axis.title.y=element_text(size=rel(1.5)))+
  theme(legend.text=element_text(size=rel(1.2)), legend.title=element_text(size=rel(1.2)))+
  scale_color_manual(values = c("#BBBBBB", "#4477AA", "#EE6677", "#228833"))
dev.off()
