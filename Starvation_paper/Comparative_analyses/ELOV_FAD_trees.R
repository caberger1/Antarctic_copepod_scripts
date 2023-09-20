library(dplyr)
library(phytools)
library(tidytree)
library(ggplot2)
library(ggtree)
library(viridis)

setwd("F:/Antarctic_copepods/trees/")

ELOV_36 <- read.tree("ELOV_36_trim.tree")
ELOV_125 <- read.tree("ELOV_125_trim.tree")
FAD <- read.tree("FAD_trim.tree")

#write.csv(ELOV_36$tip.label, file="ELOV_36_annots.csv", row.names = F, quote = F)
#write.csv(ELOV_125$tip.label, file="ELOV_125_annots.csv", row.names = F, quote = F)
#write.csv(FAD$tip.label, file="FAD_annots.csv", row.names = F, quote = F)

Gene_info <-  read.csv("ELOV_36_annots.csv")
ELOV_36 <- full_join(ELOV_36, Gene_info, by="label")
Gene_info <-  read.csv("ELOV_125_annots.csv")
ELOV_125 <- full_join(ELOV_125, Gene_info, by="label")
Gene_info <-  read.csv("FAD_annots.csv")
FAD <- full_join(FAD, Gene_info, by="label")


#propagate to internal nodes
label_internal_nodes <- function(tree){
  tree <- as_tibble(tree)
  
  internal_nodes = unique(tree$parent)
  
  #return TRUE if all descendant nodes are in the foreground
  is_foreground <- function(node){
    desc = offspring(tree, node)
    
    return( all(desc$Foreground, na.rm=T) )
  }
  
  tree$Foreground[internal_nodes] <- sapply(internal_nodes, is_foreground)
  return(tree)
}

ELOV_36 <- label_internal_nodes(ELOV_36)
ELOV_36 <- as.treedata(ELOV_36)

ELOV_125 <- label_internal_nodes(ELOV_125)
ELOV_125 <- as.treedata(ELOV_125)

FAD <- label_internal_nodes(FAD)
FAD <- as.treedata(FAD)

p <- ggtree(ELOV_36, aes(color=Foreground)) 

p + geom_text(aes(label=node), hjust=-.3) #see internal node labels

colors <- c("darkgrey", "#381A61FF", "#88A0DCFF", "#7C4B73FF")

p + geom_tiplab(aes(label=Species), color="black") + geom_tippoint(size=4, aes(color=DE)) +
  scale_color_manual(values=c("black", "#00000000", "#619CFF", "red")) 

p + geom_tiplab(aes(label=Species), color="black") + geom_tippoint(size=4, shape=21, aes(color=DE, fill=LFC)) +
  scale_color_manual(values=c("black", "#00000000", "#619CFF", "red")) + 
  scale_fill_viridis()
  
p <- ggtree(ELOV_125, aes(color=Foreground)) 

p + geom_tiplab(aes(label=Species), color="black") + geom_tippoint(size=4, aes(color=DE)) +
  scale_color_manual(values=c("black", "#00000000", "#619CFF", "red")) 



##Make FAD figure(s)
p <- ggtree(FAD, aes(color=Foreground)) 

p <- p + geom_tiplab(aes(label=node), color="black", align=TRUE, size=2) + geom_tippoint(size=4, aes(color=DE)) +
  scale_color_manual(values=c("black", "#00000000", "#619CFF", "red")) 
p

viewClade(p, MRCA(p, 71, 84))
viewClade(p, MRCA(p, 142, 84))

p <- ggtree(FAD, aes(color=Foreground)) 
p <- p + geom_tiplab(aes(label=Species), color="black", align=TRUE, size=2) + #geom_tippoint(size=4, aes(color=DE)) +
  scale_color_manual(values=c("black",  "#619CFF", "red")) 
viewClade(p, MRCA(p, 142, 84))

p <- p + geom_hilight(mapping=aes(subset = Species %in% c("Cacu", "Cpro"), fill=Species), alpha = .8, extend=0.2) +
  scale_fill_manual(values=c("#A9D18E", "#8FAADC"), labels=c("Calanoides acutus", "Calanus propinquus"))

png(res=400, file="FAD_LEGEND.png", units='in', width=10, height=8)
p + theme(text=element_text(size=20))
dev.off()

viewClade(p, MRCA(p, 142, 84))

png(res=400, file="FAD_test.png", units='in', width=4, height=12)
p <- ggtree(FAD, aes(color=Foreground)) 
p <- p + geom_tiplab(aes(label=Species), color="black", align=TRUE, size=2) +
  scale_color_manual(values=c("black", "#619CFF")) 
p <- p + theme(legend.position="none")
viewClade(p, 157, xmax_adjust = 0.1)
dev.off()

png(res=400, file="FAD_Zoom1_v2.png", units='in', width=4, height=8)
p <- ggtree(FAD, aes(color=Foreground)) 
p <- p + geom_tiplab(aes(label=Species), color="black", align=TRUE, size=2.5) +
  scale_color_manual(values=c("black", "#619CFF")) 
p <- p + theme(legend.position="none")
viewClade(p, MRCA(p, 71, 84), xmax_adjust = 0.1)
dev.off()

png(res=400, file="FAD_Zoom2.png", units='in', width=4, height=12)
p <- ggtree(FAD, aes(color=Foreground)) 
p <- p + geom_tiplab(aes(label=Species), color="black", align=TRUE, size=3) +
  scale_color_manual(values=c("black", "#619CFF")) 
p <- p + theme(legend.position="none")
viewClade(p, MRCA(p, 142, 84), xmax_adjust = 0.1)
dev.off()

##Make ELOV figure(s)
p <- ggtree(ELOV_36, aes(color=Foreground)) 

p <- p + geom_tiplab(aes(label=node), color="black", align=TRUE, size=2) + geom_tippoint(size=4, aes(color=DE)) +
  scale_color_manual(values=c("black", "#00000000", "#619CFF", "red")) 
p

p <- ggtree(ELOV_36, aes(color=Foreground)) 
p <- p + geom_tiplab(aes(label=Species), color="black", align=TRUE, size=2) + geom_tippoint(size=4, aes(color=DE)) +
  scale_color_manual(values=c("black", "#00000000", "#619CFF", "red")) 

p <- p + geom_hilight(mapping=aes(subset = Species %in% c("Cacu", "Cpro"), fill=Species), alpha = .8, extend=0.2) +
  scale_fill_manual(values=colors[3:4])

png(res=400, file="ELOV_36_v2.png", units='in', width=4, height=8)
p <- ggtree(ELOV_36, aes(color=Foreground)) 
p <- p + geom_tiplab(aes(label=Species), color="black", align=TRUE, size=2.5) +
  scale_color_manual(values=c("black", "#619CFF")) 
p <- p + theme(legend.position="none")
viewClade(p, 63, xmax_adjust = 0.1)
dev.off()

### ELOV125
p <- ggtree(ELOV_125, aes(color=Foreground)) 

p <- p + geom_tiplab(aes(label=node), color="black", align=TRUE, size=2) + geom_tippoint(size=4, aes(color=DE)) +
  scale_color_manual(values=c("black", "#00000000", "#619CFF", "red")) 
p

viewClade(p, MRCA(p, 186, 488), xmax_adjust = 0.1)
#MRCA node is 250
viewClade(p, MRCA(p, 186, 92), xmax_adjust = 0.1)
#Clade of interest is 315

p <- ggtree(ELOV_125, aes(color=Foreground)) 
p <- p + geom_tiplab(aes(label=Species), color="black", align=TRUE, size=2) + geom_tippoint(size=4, aes(color=DE)) +
  scale_color_manual(values=c("black", "#00000000", "#619CFF", "red")) 

p <- p + geom_hilight(mapping=aes(subset = Species %in% c("Cacu", "Cpro"), fill=Species), alpha = .8, extend=0.2) +
  scale_fill_manual(values=colors[3:4])
viewClade(p, MRCA(p, 186, 92), xmax_adjust = 0.1)

png(res=400, file="ELOV_125_test.png", units='in', width=4, height=12)
p <- ggtree(ELOV_125, aes(color=Foreground)) 
p <- p + geom_tiplab(aes(label=Species), color="black", align=TRUE, size=1.5) +
  scale_color_manual(values=c("black", "#619CFF")) 
p <- p + theme(legend.position="none")
viewClade(p, 250, xmax_adjust = 0.1)
dev.off()

png(res=400, file="ELOV_125_Zoom.png", units='in', width=4, height=12)
p <- ggtree(ELOV_125, aes(color=Foreground)) 
p <- p + geom_tiplab(aes(label=Species), color="black", align=TRUE, size=2.25) +
  scale_color_manual(values=c("black", "#619CFF", "#00000000", "")) 
p <- p + theme(legend.position="none")
viewClade(p, 315, xmax_adjust = 0.1)
dev.off()


##Species tree
Species <- read.tree("C20_T500.treefile")
Species <- midpoint.root(Species)

Sp_info <-  read.csv("Species_annots.csv")
Species <- full_join(Species, Sp_info, by="label")

Species <- label_internal_nodes(Species)
Species <- as.treedata(Species)

focal = c("Calanoides acutus", "Calanus propinquus")

png(res=400, units='in', height=8, width=6, file="Species_tree_Partition.png")
p <- ggtree(Species, aes(color=Foreground), size=1) + geom_tiplab(aes(label=Species_long), color="black")
p <- p + theme(legend.position = "none") + xlim(c(0,2.5)) + scale_color_manual(values=c("black", "#619CFF"))
p + geom_tippoint(size=3, shape=21, aes(fill=Group)) + scale_fill_manual(values=c("#00000000", "#8FAADC", "#A9D18E")) 
dev.off()

png(res=400, units='in', height=8, width=6, file="Species_tree2_Partition.png")
p <- ggtree(Species, aes(color=Foreground), size=1) + geom_tiplab(aes(subset=Species_long %in% focal,
  label=Species_long, fill=Group), color="black", geom="label",
                                                                  label.padding = unit(0.07, "lines"),
                                                                  label.size = 0,
  fontface="bold")                                                      
p <- p + geom_tiplab(aes(subset=!(Species_long %in% focal),
                         label=Species_long, fill=Group), color="black", geom="label",
                     label.padding = unit(0.1, "lines"),
                     label.size = 0) 
  p <- p + theme(legend.position = "none") + xlim(c(0,2.5)) + scale_color_manual(values=c("black", "#619CFF"))
p + scale_fill_manual(values=c("#00000000", "#8FAADC", "#A9D18E")) 
dev.off()

p + geom_nodelab()

p + geom_tiplab(aes(label=Species_long), color="black") + geom_tippoint(size=2, aes(fill=Group)) +
  scale_color_manual(values=c("black", "#619CFF", "#00000000", "#A9D18E", "#8FAADC")) 

p <- p + geom_tiplab(aes(label=Species_long), color="black") + geom_tippoint(size=4, aes(color=Group)) +
  scale_color_manual(values=c("#00000000", "black", "#8FAADC", "#619CFF", "#A9D18E")) 


  p <- p + geom_hilight(mapping=aes(subset = Species_long %in% c("Calanoides acutus", "Calanus propinquus"), fill=Species_long), alpha = .8, extend=0.2) +
  scale_fill_manual(values=c("#A9D18E", "#8FAADC"))
