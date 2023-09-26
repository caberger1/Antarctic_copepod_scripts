##Script tests for associations of WGCNA modules (Calanus propinquus) with environmental measurements
#throughout script replace "yourpath" with path to file

library(WGCNA)
library(psych)
library(multcomp)
library(limma)

rm(list = ls())
disableWGCNAThreads()
setwd("yourpath")

counts = read.table("limma_CPM.txt", check.names = FALSE)
WGCNA_counts <- t(counts)
options(stringsAsFactors = FALSE);
lnames = load(file = "prop_signed_bicor.RData");

n_genes = ncol(WGCNA_counts)
n_samples = nrow(WGCNA_counts)
MEs0 = moduleEigengenes(WGCNA_counts, module_colors)$eigengenes
MEs = orderMEs(MEs0)

Design <- read.csv("Propinquus_factors.csv", header = TRUE, stringsAsFactors = TRUE)
Design$Station <- factor(Design$Station)
MEs_field <- MEs[c(1:16),]

#Intercept is expression at the "reference site", in this case S1. 
options(na.action = 'na.pass');
mm = model.matrix( ~1 + Site, data = Design)

model.lm <- lapply(MEs_field, function(x){
  lm(x ~ mm)
})

M <- mm[,-1]
model.lm <- lapply(model.lm, function(x){
  update(x, ~M)
})

colnames(mm)[1] <- "Intercept"

Site_contrast <- makeContrasts( #get difference between group and mean of other 3 sites (note that the intercept is within each of the coefficients, hence the algebra)
  S1 = -(SiteS2 + SiteS3 + SiteS4)/3, 
  S2 = (SiteS2 - (SiteS3 + SiteS4)/3),  
  S3 = (SiteS3 - (SiteS2 + SiteS4)/3),  
  S4 = (SiteS4 - (SiteS2 + SiteS3)/3),  
  levels = colnames(mm))

my_contrast_fun <- function(x){
  tests = sapply(1:ncol(Site_contrast), function(i){
    t <- glht(x, linfct = t(Site_contrast[,i]))
    
    return(summary(t)$test)
  })
  
  colnames(tests) = c("Site1", "Site2", "Site3", "Site4")
  
  coef = as.data.frame(t(tests)[,c(3,6)])
  coef$padj = p.adjust(coef$pvalues, method="BH")
  
  return(coef)
}

model_contrasts = lapply(model.lm, my_contrast_fun)

mat <- lapply(model_contrasts, function(x){
  text_vect <- x$coefficients
  return(as.numeric(text_vect))
})

mat <- do.call(rbind, mat)

colnames(mat) = c("000.100", "200.000", "100.180", "-100.100")
mat <- mat[, c(2, 3, 1, 4)]


p_values = lapply(model_contrasts, function(x){
  return(x$pvalues)
})
p_values <- do.call(rbind, p_values)
p_values <- p_values[, c(2, 3, 1, 4)]

##
##remove grey module (which is not a real module) and correct for multiple testing within each contrast
p_values = p_values[-nrow(p_values),]
p_adj <- apply(p_values, 2, p.adjust, method = "BH")

mat = mat[-nrow(mat), ]

#text_mat = matrix(paste(signif(as.numeric(mat), 2), ", p=", signif(p_adj, 2), sep = ""), nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))

text_mat = matrix(nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
text_mat[which(p_adj < 0.05)] = paste(signif(as.numeric(mat), 2), ", p=", signif(p_adj, 2), sep = "")


text_mat[which(p_adj < 0.05)] = paste(text_mat[which(p_adj < 0.05)], "*")
text_mat[which(p_adj < 0.01)] = paste(text_mat[which(p_adj < 0.01)], "*", sep='')
text_mat[which(p_adj < 0.001)] = paste(text_mat[which(p_adj < 0.001)], "*", sep='')

mod_names <- substring(rownames(text_mat), 3)

#export figure as desired
#pdf(width=12, height=12, file="filename.pdf")
#png(units='in', width=12, height=12, res=600, file="filename.png")
par(mai = c(0.5, 2.5, 0.5, 0.5));

labeledHeatmap(Matrix = mat,
               xLabels = colnames(text_mat),
               yLabels = rownames(text_mat),
               ySymbols = mod_names,
               xLabelsAngle = 0,
               colorLabels = TRUE,
               setStdMargins = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_mat,
               cex.text = 0.7,
               zlim = c(-0.7,0.7),
               main = "Module-treatment relationships"
)
#dev.off()






