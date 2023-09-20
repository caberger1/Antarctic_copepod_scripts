setwd("F:/Antarctic_copepods/Acutus/limma")

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)
library(tximport)
library(readr)
library(edgeR)
library(limma)
library(pheatmap)
library(stringr)
library(adegenet)

counts_dir = "F:/Antarctic_copepods/Acutus/counts/"
count_files <- file.path(counts_dir, list.files(counts_dir))
names(count_files) <- basename(gsub(count_files, pattern = "_quant.sf", replacement = ''))

tx2gene = readr::read_tsv("F:/Antarctic_copepods/Acutus/acutus_clusters.txt")

#design
group_table = read.table("F:/Antarctic_copepods/Acutus/group_table.txt", header=TRUE)

#remove outlier (Fig. S1)
group_table = as.data.frame(group_table[-40,])
names(group_table) = "group"
count_files = count_files[-40]

txi <- tximport(count_files, type="salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
rownames(group_table) <- colnames(txi$counts)
y <- DGEList(txi$counts)

keep <- rowSums( y$counts >= 15 ) >= 4
#keep 59,540

y = y[keep,, keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
cpms = cpm(y, log=TRUE)
write.table(cpms, "Acutus_limma_CPM.txt", quote = FALSE, sep = "\t")
cpms <- read.table("Acutus_limma_CPM.txt")

field_cols <- c("#AA3377", "#66CCEE", "#CCBB44", "#228833", "#4477AA", "#5d5d5d")
field_labs <- paste0("Field, ", c("616.040", "400.040", "200.040", "200.000", "000.100", "-100.100"))
field_breaks <- c("S117", "S27", "S45", "S62", "S77", "S81")
starve_cols <- c("#E76254FF", "#EF8A47FF", "#AADCE0FF", "#72BCD5FF", "#376795FF")
starve_breaks <- c("S9", "S4", "R9", "F4", "F9")
starve_labs <- c("Starved D9", "Starved D5", "Refed D9", "Fed D5", "Fed D9")
Shape_labs <- c( rep(3, 15), rep(16, 8), rep(3, 15), rep(16, 11))


#to match Ann's figure: 
png(res=600, width=6, height=4, units='in', file="F:/Antarctic_copepods/starve/figs/Acutus_PCA_full.png")
pca <- prcomp(t(cpms))
df <- cbind(group_table, pca$x)
PC1_var <- summary(pca)$importance[2]
PC2_var <- summary(pca)$importance[5]
p <- ggplot(df) + geom_point(aes(x=PC1, y=PC2, color=group, shape=group), size=3)+ 
  ggtitle("Acutus PCA, CPM") + xlab(paste0("PC1 (", signif(PC1_var*100, 3), "%)")) + 
  ylab(paste0("PC2 (", signif(PC2_var*100, 3), "%)")) + theme_classic()

p + scale_color_manual(name = "Group",
                       values = c(field_cols, starve_cols), 
                       breaks = c(field_breaks, starve_breaks),
                       labels = c(field_labs, starve_labs)) + 
  scale_shape_manual(name = "Group", labels = c(field_labs, starve_labs),
                     breaks = c(field_breaks, starve_breaks), 
                     values = c(rep(3, 6), rep(16, 5))) +
  scale_y_reverse() + theme(text = element_text(size = 20), plot.title = element_blank())
dev.off()

#And just the Expt
png(res=600, width=6, height=4, units='in', file="F:/Antarctic_copepods/starve/figs/Acutus_PCA_Expt.png")
pca <- prcomp(t(cpms)[-c(1:15, 24:38),])
df <- cbind(group_table[-c(1:15, 24:38),,drop=F], pca$x)
PC1_var <- summary(pca)$importance[2]
PC2_var <- summary(pca)$importance[5]
p <- ggplot(df) + geom_point(aes(x=PC1, y=PC2, color=group), size=3)+ 
  ggtitle("Acutus PCA, CPM") + xlab(paste0("PC1 (", signif(PC1_var*100, 3), "%)")) + 
  ylab(paste0("PC2 (", signif(PC2_var*100, 3), "%)")) + theme_classic()

p + scale_color_manual(name = "Group",
                       values = starve_cols, 
                       breaks = starve_breaks,
                       labels = starve_labs) +
  scale_x_reverse() + scale_y_reverse() + theme(text = element_text(size = 20), plot.title = element_blank())
dev.off()

#design
treatments <- read.table("../sample_table.txt", header=TRUE, stringsAsFactors = TRUE)
treatments$site = factor(treatments$site, levels=c("45" ,"27"  ,"62"  ,"77" ,"81" , "117"))
treatments$day = factor(treatments$day)
treatments = treatments[-40,]

options(na.action = 'na.pass');
mm = model.matrix( ~0+ group + food*day + site, data = treatments)
mm[is.na(mm)] <- 0

design <- mm[,colSums(mm) != 0]

##create contrast matrix
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- gsub("food", "", colnames(design))
colnames(design) <- gsub(":", "_", colnames(design))

#we can use DAPC to show that "refed" reverses starved. 
#code modified from Groves et al. 2015
fed_dat = cpms[,c(16:23,39:41,46:49)]
refed_dat = cpms[,c(42:45)]

dp <- dapc(t(fed_dat), treatments[!is.na(treatments$food),]$food[-c(12:15)], n.da=1, perc.pca=80)

assemble_dp_coords = function(dp.object){
  res = tibble(Run = rownames(dp.object$ind.coord),
               LD1 = dp.object$ind.coord[,'LD1'],
               group = as.character(dp.object$grp))
}

coords = assemble_dp_coords(dp)
coords %>% 
  ggplot(aes(x=LD1,fill=group)) +
  geom_density(alpha=0.8) + theme_bw() + xlim(c(-10, 10)) + geom_point( aes(y = 0))

assemble_pred_coords = function(pred.object, group.vector){
  res = tibble(Run = rownames(pred.object$ind.scores),
               LD1 = pred.object$ind.scores[,'LD1'],
               group = group.vector)
}

pred_Refed <- predict.dapc(dp, newdata = (t(refed_dat)))
Refed_coords <- assemble_pred_coords(pred_Refed, treatments$food[c(12:15)])

df <- rbind(coords, Refed_coords)
df[is.na(df$group),]$group = "refed"

png(res=400, file="F:/Antarctic_copepods/starve/figs/DAPC_refed.png",units='in', height=4, width=4)
df %>% 
  ggplot(aes(x=LD1,fill=group)) +
  geom_density(alpha=0.7) + geom_point( aes(y = 0), color = "black", shape=21) + 
  scale_fill_manual(values=c(starve_cols[5], starve_cols[3], starve_cols[1]), labels=c("Fed", "Refed", "Starved")) + 
  theme_bw() + xlim(c(-10, 10)) + 
  theme(legend.position="top", text = element_text(size = 16), plot.title = element_blank(), legend.title = element_blank()) + ylab("Density") 
dev.off()
###

#This bit may look confusing, but I have double checked that these contrasts are correct.
#It's just the necessary algebra to deal with the uneven sample design.

cmat <- makeContrasts(
  ship_v_S1 = ( (5*lab + 3*day9) + (2*starved + starved_day9) + refed)/5 - field,
  fed_v_S1 = (2*lab + day9)/2 - field,
  ship_v_Field = ( (5*lab + 3*day9) + (2*starved + starved_day9) + refed)/5 - 
    (6*field + site27+site62+site77+site81+site117)/6,
  starved_total = (2*starved + starved_day9)/2, 
  starved_late = starved + starved_day9,
  starved_early = starved,
  starve_day_int = starved_day9,
  refed_v_fed = refed,
  refed_v_starved = refed - (starved+starved_day9),
  levels = colnames(design))

#ref level should be Site 45
vq <- voomWithQualityWeights(y, design, plot=TRUE)

vfit <- lmFit(vq, design)
vfit <- contrasts.fit(vfit, contrasts=cmat)

efit <- eBayes(vfit, robust=TRUE)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))

#Will want the standard errors and degrees of freedom for comparative analyses
write.csv(efit$stdev.unscaled*sqrt(efit$s2.post), file="F:/Antarctic_copepods/starve/acutus/SE_acutus.csv")
write.csv(cbind(rownames(efit$coefficients), sqrt(efit$s2.post)), file="F:/Antarctic_copepods/starve/acutus/SD_acutus.csv")
write.csv(cbind(rownames(efit$coefficients), efit$df.total), file="F:/Antarctic_copepods/starve/acutus/df_acutus.csv")

lapply(seq(1,8), function(i){
  write.csv(topTable(efit, coef=i, confint=TRUE, number = Inf), file=paste0(colnames(efit)[i], "_efit.csv"))
})

##WGCNA
library(WGCNA)
library(psych)
library(multcomp)

disableWGCNAThreads()
load("F:/Antarctic_copepods/Acutus/WGCNA/acutus_signed_bicor.RData")

#filter low-expression
counts = read.table("Acutus_limma_CPM.txt", check.names = FALSE)

##WGCNA wants genes as columns
dat_expr <- t(counts)

#checks for missing entries and zero-variance genes
gsg = goodSamplesGenes(dat_expr, verbose = 3);
gsg$allOK

##Cluster samples to look for outliers
sample_tree = hclust(dist(dat_expr), method = "average");
plot(sample_tree, main = "Sample clustering to detect outliers", sub="", xlab='')

## Choose soft thresholding exponent (for calculating adjacency)

powers = c(c(1:10), seq(from=12, to = 20, by=2))

sft <- pickSoftThreshold(dat_expr, powerVector=powers, verbose=5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red");

abline(h=0.90,col="red")

soft_power = 7 

options(stringsAsFactors = FALSE);
lnames = load(file = "F:/Antarctic_copepods/Acutus/WGCNA/acutus_signed_bicor.RData");

## want to identify modules associated with treatment effects
n_genes = ncol(dat_expr)
n_samples = nrow(dat_expr)

MEs0 = moduleEigengenes(dat_expr, module_colors)$eigengenes
MEs = orderMEs(MEs0)

treatments <- read.table("F:/Antarctic_copepods/Acutus/sample_table.txt", header = TRUE, stringsAsFactors = TRUE)

treatments$site = factor(treatments$site, levels=c("45" ,"27"  ,"62"  ,"77" ,"81" , "117"))
treatments$day = factor(treatments$day)
treatments = treatments[-40,]


options(na.action = 'na.pass');
mm = model.matrix( ~ group + food*day + site, data = treatments)
mm[is.na(mm)] <- 0


##categorical data
model.lm <- lapply(MEs, function(x){
  lm(x ~ mm, drop.unused.levels = TRUE)
})

##the problem here is that we want to make specific comparisons that are not in this table of contrasts. as it is now, 
##each factor has a reference level. Meaning, "Lab" is actually the difference betweeen lab samples and Site27, where we really want to average
##over all field sites. Similarly, "refed" and "starved" both compare versus "fed" ON DAY 4. We can do this with a generalied linear hypothesis test,
#using a vector of column coefficients. Let's consult cmat from the DE analysis, which basically already did this.

##we need to remove the NA model coefficients. These are the intercept and refed:day9, so

attr.assign = attr(mm, "assign")
attr.contrasts = attr(mm, "contrasts")
##these are the NA column numbers
rmvCols = c(1, 11)
M = mm[,-rmvCols] # Remove intercept column (required) and NA columns
# Reassign the "assign" attribute, removing the corresponding elements from it.
attr(M, "assign") = attr.assign[-rmvCols]
# Reassign the "contrasts" attribute to its original value
attr(M, "contrasts") = attr.contrasts

model.lm <- lapply(model.lm, function(x){
  update(x, ~M)
})

design <- mm[,colSums(mm) != 0]

##create contrast matrix
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- gsub("food", "", colnames(design))
colnames(design) <- gsub(":", "_", colnames(design))

##build up my contrasts manuallly...
cmat <- makeContrasts(
  lab_v_S1 = ( (5*lab + 3*day9) + (2*starved + starved_day9) + refed)/5,
  fed_v_S1 = (2*lab + day9)/2,
  starved_total = (2*starved + starved_day9)/2, 
  starved_early = starved,
  starved_late = starved + starved_day9,
  starve_day_int = starved_day9,
  refed_v_fed = refed,
  site27 = site27,
  site62 = site62, 
  site77 = site77,
  site81 = site81,
  site177 = site117,
  levels = colnames(design))

#drop refed_day9
cmat = cmat[-11,]

my_contrast_fun <- function(x){
  
  tests = sapply(1:ncol(cmat), function(i){
    t <- glht(x, linfct = t(cmat[,i]))
    
    return(summary(t)$test)
  })
  
  colnames(tests) = c("Ship", "Fed_v_S1", "Starve", "Early", "Late", "Int", "Refed", "Site 27", "Site 62", "Site 77", "Site 81", "Site 177")
  
  coef = as.data.frame(t(tests)[,c(3,6)])
  coef$padj = p.adjust(coef$pvalues, method="BH")
  
  return(coef)
}

test = glht(model.lm[[1]], linfct = t(cmat[,2]))

model_contrasts = lapply(model.lm, my_contrast_fun)

##heatmap of fitted values ("means"), and below each row the Pearson coefficient

mat <- lapply(model_contrasts, function(x){
  text_vect <- x$coefficients
  #text_vect <- text_vect[c(1:11)]
  return(as.numeric(text_vect))
})

mat <- do.call(rbind, mat)
colnames(mat) = c("Ship", "Fed_v_S1", "Starve", "Early", "Late", "Int", "Refed", "Site 27", "Site 62", "Site 77", "Site 81", "Site 177")

p_values = lapply(model_contrasts, function(x){
  return(x$pvalues)
})
p_values <- do.call(rbind, p_values)

##remove grey module (which is not a real module) and correct for multiple testing within each variable (model term)
p_values = p_values[-nrow(p_values),]
p_adj <- apply(p_values, 2, p.adjust, method = "BH")

mat = mat[-nrow(mat), ]

text_mat = matrix(paste(signif(as.numeric(mat), 2), ", p=", signif(p_adj, 2), sep = ""), nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))

text_mat[which(p_adj < 0.05)] = paste(text_mat[which(p_adj < 0.05)], "*")
text_mat[which(p_adj < 0.01)] = paste(text_mat[which(p_adj < 0.01)], "*", sep='')
text_mat[which(p_adj < 0.001)] = paste(text_mat[which(p_adj < 0.001)], "*", sep='')

mod_names <- substring(rownames(text_mat), 3)

##plot for Supplement (relabel columns to match manuscript table)
mat <- mat[,c(3,4,5,6,7,2)]
text_mat <- text_mat[,c(3,4,5,6,7,2)]
colnames(text_mat) <- c("Starved (overall)", "Starved (Day 5)", "Starved (Day 9)", 
                        "Starve:Day interaction", "Refed vs. Fed", "Expt. vs. Field")

sizeGrWindow(20,12)

png(units='in', width=20, height=12, res=600, file="Acutus_WGCNA_Supp.png")
par(mai = c(0.5, 2, 0.5, 0.2));

labeledHeatmap(Matrix = mat,
               xLabels = colnames(text_mat),
               yLabels = rownames(text_mat),
               ySymbols = mod_names,
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               xColorOffset = 0.01,
               colorLabels = TRUE,
               setStdMargins = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_mat,
               cex.text = 1.25,
               cex.lab.x = 1.75,
               cex.lab.y = 1.25,
               zlim = c(-0.5,0.5),
               
               main = "Module-treatment relationships"
)
dev.off()
