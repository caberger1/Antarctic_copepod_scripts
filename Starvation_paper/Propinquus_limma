setwd("F:/Antarctic_copepods/Propinquus/salmon/")

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
library(paletteer)

#read in count files
counts_dir = "F:/Antarctic_copepods/Propinquus/salmon/counts/"
count_files <- file.path(counts_dir, list.files(counts_dir))
names(count_files) <- basename(gsub(count_files, pattern = "_quant.sf", replacement = ''))

tx2gene = readr::read_tsv("F:/Antarctic_copepods/Propinquus/prop_clusters.txt")

treatments <- read.table("F:/Antarctic_copepods/Propinquus/salmon/sample_table.txt", header=TRUE)
treatments$site = factor(treatments$site)
treatments$day = factor(treatments$day)

#We remove an outlier (Fig. S2)
treatments = treatments[-24,]
count_files = count_files[-24]

txi <- tximport(count_files, type="salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
TPM10k <- (txi$abundance * nrow(txi$abundance) )/ 10**4
write.csv(file="Prop_TPM10k.txt")

y <- DGEList(txi$counts)
keep <- rowSums( y$counts >= 15 ) >= 4
#retains 76,556 genes

y = y[keep,, keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
cpms = cpm(y, log=TRUE)
write.table(cpms, "prop_limma_CPM.txt", quote = FALSE, sep = "\t")
cpms <- read.table("prop_limma_CPM.txt")

field_cols <- c("#4477AA", "#228833", "#EE6677", '#5d5d5d')
field_labs <- paste0("Field, ", c("000.100", "200.000", "100.180", "-100.100"))
field_breaks <- c("S1", "S2", "S3", "S4")
starve_cols <- c("#E76254FF", "#EF8A47FF", "#72BCD5FF", "#376795FF")
starve_breaks <- c("S9", "S5", "F5", "F9")
starve_labs <- c("Starved D9", "Starved D5", "Fed D5", "Fed D9")

group_table <- data.frame(group = c(rep(c("S1", "S2", "S3", "S4"), each=4), 
                                       rep("F5", 4), rep("S5", 3), rep(c("F9", "S9"), each=4)))
#to match Ann's figure: 
png(res=600, width=6, height=4, units='in', file="F:/Antarctic_copepods/starve/figs/Prop_PCA_full_080123.png")
pca <- prcomp(t(cpms))
df <- cbind(group_table, pca$x)
PC1_var <- summary(pca)$importance[2]
PC2_var <- summary(pca)$importance[5]
p <- ggplot(df) + geom_point(aes(x=PC1, y=PC2, color=group, shape=group), size=3)+ 
  ggtitle("Prop PCA, CPM") + xlab(paste0("PC1 (", signif(PC1_var*100, 3), "%)")) + 
  ylab(paste0("PC2 (", signif(PC2_var*100, 3), "%)")) + theme_classic()  + theme(text = element_text(size = 20), plot.title = element_blank())

p + scale_color_manual(name = "Group",
                       values = c(field_cols, starve_cols), 
                       breaks = c(field_breaks, starve_breaks),
                       labels = c(field_labs, starve_labs)) + 
  scale_shape_manual(name = "Group", labels = c(field_labs, starve_labs),
                     breaks = c(field_breaks, starve_breaks), 
                     values = c(rep(3, 4), rep(16, 4))) 
dev.off()

#And just the Expt
png(res=600, width=6, height=4, units='in', file="F:/Antarctic_copepods/starve/figs/Prop_PCA_Expt.png")
pca <- prcomp(t(cpms)[-c(1:16),])
df <- cbind(group_table[-c(1:16),,drop=F], pca$x)
PC1_var <- summary(pca)$importance[2]
PC2_var <- summary(pca)$importance[5]
p <- ggplot(df) + geom_point(aes(x=PC1, y=PC2, color=group), size=3)+ 
  ggtitle("Prop PCA, CPM") + xlab(paste0("PC1 (", signif(PC1_var*100, 3), "%)")) + 
  ylab(paste0("PC2 (", signif(PC2_var*100, 3), "%)")) + theme_classic()

p + scale_color_manual(name = "Group",
                       values = starve_cols, 
                       breaks = starve_breaks,
                       labels = starve_labs) +
  scale_y_reverse() + scale_x_reverse()  + theme(text = element_text(size = 20), plot.title = element_blank())
dev.off()



treatments$site <- relevel(treatments$site, ref="S1")


options(na.action = 'na.pass');
mm = model.matrix( ~0+ group + food*day + site, data = treatments)
mm[is.na(mm)] <- 0

design <- mm[,colSums(mm) != 0]

##create contrast matrix
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- gsub("food", "", colnames(design))
colnames(design) <- gsub(":", "_", colnames(design))

# We use dummy coding to create the design matrix, and then have to use some algebra to get the contrasts of interest.
# The interpretation of the coefficients in the design column is: lab is the mean at F4. Field is the mean at S1 (reference levels).
# Starved coefficient is the difference between starve and fed ON DAY 4 (so "early"). Starve_day9 is the additional effect of day9 on the starve effect
# To get the overall starve effect, we have to average these. A similar logic gives the average of Fed4 + Fed9, and then we add half of "starved_total" to that to end up with the average of all the lab/ship samples.

cmat <- makeContrasts(
  lab_v_S1 = (2*lab + day9)/2 + (2*starved + starved_day9)/4 - field,
  fed_v_S1 = (2*lab + day9)/2 - field,
  starved_total = (2*starved + starved_day9)/2, 
  starved_late = starved + starved_day9,
  starved_early = starved,
  starve_day_int = starved_day9,
  levels = colnames(design))

vq <- voomWithQualityWeights(y, design, plot=TRUE)

vfit <- lmFit(vq, design)
vfit <- contrasts.fit(vfit, contrasts=cmat)

efit <- eBayes(vfit, robust=TRUE)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))

#Will want the standard errors, SD, and degrees of freedom for comparative analyses
write.csv(efit$stdev.unscaled*sqrt(efit$s2.post), file="F:/Antarctic_copepods/starve/propinquus/SE_prop_080123.csv")
write.csv(cbind(rownames(efit$coefficients), sqrt(efit$s2.post)), file="F:/Antarctic_copepods/starve/propinquus/SD_prop_080123.csv")
write.csv(cbind(rownames(efit$coefficients), efit$df.total), file="F:/Antarctic_copepods/starve/propinquus/df_prop_080123.csv")

lapply(seq(1,6), function(i){
  write.csv(topTable(efit, coef=i, confint=TRUE, number = Inf), file=paste0(colnames(efit)[i], "_efit_080123.csv"))
})

##WGCNA
library(WGCNA)
library(psych)
disableWGCNAThreads()

counts = read.table("prop_limma_CPM.txt", check.names = FALSE)

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

abline(h=0.85,col="red")

soft_power = 7

options(stringsAsFactors = FALSE);
lnames = load(file = "F:/Antarctic_copepods/Propinquus/prop_signed_bicor.RData");

## want to identify modules associated with treatment effects
n_genes = ncol(dat_expr)
n_samples = nrow(dat_expr)

MEs0 = moduleEigengenes(dat_expr, module_colors)$eigengenes
MEs = orderMEs(MEs0)

treatments <- read.table("F:/Antarctic_copepods/Propinquus/WGCNA/CTD_table.txt", header = TRUE)
treatments$group = factor(treatments$group)
treatments$food = factor(treatments$food)
treatments$site = factor(treatments$site)
treatments$day = factor(treatments$day)


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

##we need to remove the NA model coefficients. These are the intercept, so

attr.assign = attr(mm, "assign")
attr.contrasts = attr(mm, "contrasts")
##these are the NA column numbers
rmvCols = c(1)
M = mm[,-rmvCols] # Remove intercept column (required) and and other NA columns
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
  lab_v_S4 = (2*lab + day9)/2 + (2*starved + starved_day9)/4,
  fed_v_S4 = (2*lab + day9)/2,
  starved_total = (2*starved + starved_day9)/2, 
  starved_early = starved,
  starved_late = starved + starved_day9,
  starve_day_int = starved_day9,
  #site_S2 = siteS2,
  #site_S3 = siteS3,
  #site_S4 = siteS4,
  levels = colnames(design))

cmat = cmat[-11,]

my_contrast_fun <- function(x){
  

  
  tests = sapply(1:ncol(cmat), function(i){
    t <- glht(x, linfct = t(cmat[,i]))
    
    return(summary(t)$test)
  })
  
  colnames(tests) = c("Ship", "Fed_v_S1", "Starve", "Early", "Late", "Int")
  
  coef = as.data.frame(t(tests)[,c(3,6)])
  coef$padj = p.adjust(coef$pvalues, method="BH")
  
  return(coef)
}

library(multcomp)

test = glht(model.lm[[1]], linfct = t(as.matrix(cmat[,2])))

model_contrasts = lapply(model.lm, my_contrast_fun)

##heatmap of fitted values ("means"), and below each row the Pearson coefficient

mat <- lapply(model_contrasts, function(x){
  text_vect <- x$coefficients
  return(as.numeric(text_vect))
})

mat <- do.call(rbind, mat)
colnames(mat) = c("Ship",  "Fed_v_S1", "Starve", "Early", "Late", "Int")


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
mat <- mat[,c(3,4,5,6,2)]
text_mat <- text_mat[,c(3,4,5,6,2)]
colnames(text_mat) <- c("Starved (overall)", "Starved (Day 5)", "Starved (Day 9)", "Starve:Day interaction", "Expt. vs. Field")

sizeGrWindow(20,12)

png(units='in', width=20, height=12, res=600, file="Prop_WGCNA_Supp.png")
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

