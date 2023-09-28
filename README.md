# Antarctic_copepod_scripts

This repo contains R code and various data files for the companion manuscripts "Comparative analysis of the molecular starvation response of Southern Ocean copepods" by Berger et al. (2023) (https://doi.org/10.1101/2023.08.17.553703) and "Nutritional condition drives spatial variation in physiology of Antarctic lipid-storing copepods" by Berger, Steinberg, & Tarrant (2023) (https://doi.org/10.1101/2023.09.25.559317). 

Within the starvation paper, there are folders for differential expression/WGCNA analyses, comparative analyses (including selection/phylogenetic analyses), and physiology (lipids and enzyme activity). 

For the Starvation experiments, the main intermediate results file is “df_060322_Annotated.csv”. This file contains detailed information for each gene in both species, including expression data, Orthogroup assignments, and BLAST annotations. “df_Ship_060322_Annotated.csv” contains the same information, except with fold-change/expression data from the Field vs. ship comparison. These files are used as input for “Comparative_Analyses.Rmd”, which contains the bulk of the comparative analyses presented in the paper. 

For variance of gene expression analysis:
“Acutus_Field_Input.csv” and “Prop_Field_Input.csv” contain raw counts, library sizes, and lengths for each gene and field sample. These files are used as input for “Overdisp_HPC_Acutus/Prop.R”, which fit negative binomial models to calculate overdispersion (as described in the Methods and Fair et al., 2020 (https://doi.org/10.7554/eLife.59929) and write the files “Acutus/Prop_Overdisp_Estimates_Field.csv”. These files are then used as input for “OverDisp_processing.R”, which calculates dispersion as the residual of a regression of overdispersion against mean expression, and writes the files “Acutus/Prop_Field_Disp.csv”. These values of dispersion are the values used in “df_060322_Annotated.csv” and analyzed within “Comparative_Analyses.Rmd”. 
