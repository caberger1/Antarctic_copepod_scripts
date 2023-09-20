##For running on HPC, calculate overdispersion etc.
library(purrr)
library(broom)
library(dplyr)
library(MASS)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]

set.seed(0)

Prop_df <- read.csv(input)

ModelExpression <- function(df) {
  glm.nb(counts ~ 1 + offset(log(LibrarySize))-offset(log(Length)), data=df)
}

#Now create a wrapper function around the modelling function to return NA in case of errors
#This allows purr::map to use my modeling function without breaking the loop when an error happens
ModelExpression.possibly <- possibly(ModelExpression, NA)

get_model_info <- function(model){
   if (is.na(model)){
	return(NA)
   }

  a <- summary(model)
  a <- data.frame(a$coefficients)
  names(a) <- c("estimate", "std.error", "statistic", "p.value")
  return( tibble(a) )  
}

ModelFits <- Prop_df %>%
  group_by(Gene) %>%
  nest() %>%
  mutate(model = map(data, ModelExpression.possibly),
         tidier = map(model, get_model_info),
         theta = map_dbl(model, "theta", .default=NA)) %>%
  unnest(c(tidier, theta))

out <- ModelFits %>% dplyr::select(Gene, estimate, std.error, statistic, p.value, theta)

write.csv(out, file="Prop_Overdisp_Estimates_Field.csv", row.names = F, quote=F)




