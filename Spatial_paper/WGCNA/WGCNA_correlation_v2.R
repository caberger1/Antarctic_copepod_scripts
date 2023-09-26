##WGCNA
library(WGCNA)
library(psych)
library(multcomp)
library(limma)
library(readr)
library(stringr)
rm(list = ls())
disableWGCNAThreads()

# I changed to read_delim to make it faster
# counts <- read.csv("/Users/atarrant/Documents/Project_Antarctic/CaLog2_for_WGCNA.csv", check.names = FALSE)
counts <- read_delim("/Users/atarrant/Documents/Project_Antarctic/CaLog2_for_WGCNA.csv", ",", col_names = TRUE)
names(counts)[names(counts) == '...1'] <- 'Samples'
counts <- counts[order(counts$Samples),] 

options(stringsAsFactors = FALSE);
lnames = load(file = "/Users/atarrant/Documents/Project_Antarctic/Ca_DE/acutus_signed_bicor.RData");

#use this to annotate transcriptome with modules
a <- as.data.frame(noquote(names(module_labels)))
a <- noquote(cbind(a, module_colors))
names(a)[names(a) == 'noquote(names(module_labels))'] <- 'cluster'
Annotation <- read.table("/Users/atarrant/Documents/Project_Antarctic/Ca_DE/Condensed_Annotation.txt", header=TRUE, fill=TRUE)
Annot_with_modules <- merge(a, Annotation, by='cluster')
#write.csv(Annot_with_modules, file="/Users/atarrant/Documents/Project_Antarctic/Ca_DE/Annot_with_modules.csv")

#Also get add WGCNA modules to annotation file of all clusters (not just the ones with Blast results)
Long_annot <- read.table("/Users/atarrant/Documents/Project_Antarctic/Ca_DE/annot_Allclusters_concise.txt", header=TRUE, fill=TRUE)
names(a)[names(a) == 'cluster'] <- 'Cluster'
All_Annot_with_modules <- merge(a, Long_annot, by='Cluster', all.y=TRUE)
#write.csv(All_Annot_with_modules, file="/Users/atarrant/Documents/Project_Antarctic/Ca_DE/annot_Allclusters_concise_with_modules.csv")
#write.table(All_Annot_with_modules, file="/Users/atarrant/Documents/Project_Antarctic/Ca_DE/annot_Allclusters_concise_with_modules.txt")


#Cory added these lines because he had filtered the genes based on field data. He had 59,540 genes. 
#First 2 lines below: We filtered out the genes (labels/colors) from his set that aren't in my set...that gave us 49005
#Third line below: I had 49026 genes...so 21 of mine weren't in his set. We filtered those out, giving us 49005
filtered_module_labels <- module_labels[names(counts)]
filtered_module_colors <- module_colors[which(names(module_labels) %in% names(counts))]

WGCNA_counts <- counts[,names(counts) %in% names(filtered_module_labels)]

# To get numbers of genes in each module also need to read in All_Annot_with_modules if you haven't run the lines
library(dplyr)
All_Annot_with_modules %>% count(module_colors) %>% 
  write.csv(file = "C:/Users/atarrant/Documents/Paper_Ant_Spatial/WGCNA/AnnotGenes_perModule.csv")

## want to identify modules associated with treatment effects
n_genes = ncol(WGCNA_counts)
n_samples = nrow(WGCNA_counts)

MEs0 = moduleEigengenes(WGCNA_counts, filtered_module_colors)$eigengenes
MEs_field = orderMEs(MEs0)


###Just do field correlations
Ac_design <- read.csv("C:/Users/atarrant/Documents/Paper_Ant_Spatial/WGCNA/Acutus_factors.csv", header = TRUE, stringsAsFactors = TRUE)
Ac_design$log_Chl <- log10(Ac_design$int_chl)
Ac_design$log_Calanoides <- log10(Ac_design$Calanoides)
Ac_design <- Ac_design[order(Ac_design$Sample),] 

#write out data for scatterplots
Sample <- counts$Samples
MEs_named <- cbind(Sample, MEs_field)
write.csv(MEs_named, file="C:/Users/atarrant/Documents/Paper_Ant_Spatial/WGCNA/Acutus_MEsField2.csv", row.names = TRUE)

colMeans(MEs_named[11:20,-1])
colMeans(MEs_named[c(1:10,21:30),-1])

module_Cor <- cor(MEs_field, Ac_design[,c(15,17,10,18)]) 
##Select the "avg_upper_temp" (10), "Site_group" (15), log_int_chl" (17),and "log_Ca_abundance" (18) columns
#ordered in a way that makes a nice figure and facilitates discussion
module_Pval <- corPvalueStudent(module_Cor, n_samples)

mat = module_Cor[-nrow(module_Cor), ] ##remove MEgrey module
p_values <- module_Pval[-nrow(module_Pval),] #remove MEgrey module
p_adj <- apply(p_values, 2, p.adjust, method = "BH")
env_p_adj <- p_adj

text_mat = matrix(paste(signif(as.numeric(mat), 2), ", p=", signif(p_adj, 2), sep = ""), nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))

text_mat[which(p_adj < 0.05)] = paste(text_mat[which(p_adj < 0.05)], "*")
text_mat[which(p_adj < 0.01)] = paste(text_mat[which(p_adj < 0.01)], "*", sep='')
text_mat[which(p_adj < 0.001)] = paste(text_mat[which(p_adj < 0.001)], "*", sep='')

mod_names <- substring(rownames(text_mat), 3)

sizeGrWindow(20,12)

png(units='in', width=20, height=12, res=600, file="C:/Users/atarrant/Documents/Paper_Ant_Spatial/WGCNA/WGCNA_Cor_Environment.png")
pdf(file="C:/Users/atarrant/Documents/Paper_Ant_Spatial/WGCNA/WGCNA_Cor_Environmental.pdf")

par(mai = c(0.5, 1.5, 0.5, 0.5));
labeledHeatmap(Matrix = mat,
               xLabels = colnames(text_mat),
               yLabels = rownames(text_mat),
               ySymbols = str_to_title(mod_names),
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               colorLabels = TRUE,
               setStdMargins = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_mat,
               cex.text = 0.7,
               zlim = c(-0.7,0.7),
               main = "Module-Factor relationships"
)
dev.off()


##and for the site-wise comparisons...

Ac_design$Site <- paste0("S", rep(1:6, each=5))
Ac_design$Site <- factor(Ac_design$Site, levels=c("S3" ,"S1"  ,"S2"  ,"S4" ,"S5" , "S6"))

#need to have a reference site, in this case S45/200.040
options(na.action = 'na.pass');
mm = model.matrix( ~1 + Site, data = Ac_design)

##categorical data
model.lm <- lapply(MEs_field, function(x){
  lm(x ~ mm)
})

#Intercept is expression at the reference site. P-value is not meaningful in that column
get_LM_output <- function(x){
  
  coef <- as.data.frame(x$coefficients[c(1,3:7)])
  coef$P_val <- summary(x)$coefficients[,4]
  names(coef)[1] <- "coefficients"
  
  return(coef)
}

model_output = lapply(model.lm, get_LM_output)

##heatmap of fitted values ("means"), and below each row the Pearson coefficient

mat <- lapply(model_output, function(x){
  text_vect <- x$coefficients
  return(as.numeric(text_vect))
})

mat <- do.call(rbind, mat)
colnames(mat) = levels(Ac_design$Site)

p_values = lapply(model_output, function(x){
  return(x$P_val)
})
p_values <- do.call(rbind, p_values)

##We can also specifically do the "Site Group" comparison:
# we want "S3" and "S4" versus everything else
#fiddle with some formatting for the intercept colum
M <- mm[,-1]
model.lm <- lapply(model.lm, function(x){
  update(x, ~M)
})

colnames(mm)[1] <- "Intercept"
Site_contrast <- makeContrasts(
  Site_groups = ( ((SiteS4 + Intercept)/2) - ((SiteS1+SiteS2+SiteS5+SiteS6)/4) ) ,
  levels=colnames(mm)
)
  
my_contrast_fun <- function(x){
  
  contrast = glht(x, linfct = t(Site_contrast[,1]))
  test <- summary(contrast)$test
  
  coef = c(test$coefficients[[1]], test$pvalues[[1]])
  
  return(coef)
}

site_contrasts <- lapply(model.lm, my_contrast_fun)

coef <- do.call(rbind,site_contrasts)[,1]
p_values <- do.call(rbind,site_contrasts)[,2]
p_values <- p_values[-length(p_values)] #dropping grey module
p_adj <- p.adjust(p_values, method="BH")




##
##remove grey module (which is not a real module) and correct for multiple testing within each variable (model term)
test <-cbind(p_adj, env_p_adj)

mat = mat[-nrow(mat), ]

text_mat = matrix(paste(signif(as.numeric(mat), 2), ", p=", signif(p_adj, 2), sep = ""), nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))

text_mat[which(p_adj < 0.05)] = paste(text_mat[which(p_adj < 0.05)], "*")
text_mat[which(p_adj < 0.01)] = paste(text_mat[which(p_adj < 0.01)], "*", sep='')
text_mat[which(p_adj < 0.001)] = paste(text_mat[which(p_adj < 0.001)], "*", sep='')

mod_names <- substring(rownames(text_mat), 3)

sizeGrWindow(20,12)
par(mai = c(0.5, 3, 0, 0));

png(units='in', width=20, height=12, res=600, file="F:/Antarctic_copepods/Acutus/Figs/WGCNA_Cor_Sites.png")

labeledHeatmap(Matrix = mat,
               xLabels = colnames(text_mat),
               yLabels = rownames(text_mat),
               ySymbols = mod_names,
               xLabelsAngle = 0,
               colorLabels = TRUE,
               setStdMargins = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_mat,
               cex.text = 0.95,
               zlim = c(-0.5,0.5),
               main = "Module-treatment relationships"
)
dev.off()