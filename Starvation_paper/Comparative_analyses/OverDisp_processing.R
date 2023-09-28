library(tximport)
library(purrr)
library(broom)
library(MASS)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(edgeR)
library(tibble)
library(FSA)
library(forcats)

Acutus_Field <- data.table( read.csv("Acutus_Overdisp_Estimates_Field.csv") )
Prop_Field <- data.table( read.csv("Prop_Overdisp_Estimates_Field.csv") )

Acutus_Field <- Acutus_Field[!is.na(theta)] #15 NA
Prop_Field <- Prop_Field[!is.na(theta)] #1 NA

Acutus_Field <- Acutus_Field[theta < 1000] #7200 blew up, n=52,325
Prop_Field <- Prop_Field[theta < 1000] #1811 blew up, n=74,744

Acutus_Field %>%
  ggplot(aes(x=estimate, y=log(1/theta))) +
  stat_binhex(bins=100) +
  xlab("mu") +
  ylab("log(1/theta)") +
  scale_fill_viridis_c(option="C") +
  geom_smooth(method="loess", method.args=list(degree=1), se=FALSE) +
  theme_bw() +
  theme(aspect.ratio = 1)

#assign the loess trend to object
Loess.Trend <- loess(log(1/theta)~estimate, data = Acutus_Field, degree = 1)

#Add a column for the residual, which we will call 'dispersion'
Acutus_Field <- Acutus_Field %>%
  ungroup() %>%
  drop_na() %>%
  mutate(Dispersion = residuals(Loess.Trend))

#Plot the residuals, with a horizontal blue line as a way to visually communicate that we are plotting the residuals from the loess trend
Acutus_Field %>%
  ggplot(aes(x=estimate, y=Dispersion)) +
  stat_binhex(bins=100) +
  xlab("mu") +
  ylab("Dispersion") +
  geom_hline(yintercept = 0, color="blue") +
  scale_fill_viridis_c(option="C") +
  theme_bw() +
  theme(aspect.ratio = 1)

write.csv(Acutus_Field, quote=F, row.names=F, file="Acutus_Field_Disp.csv")

Prop_Field %>%
  ggplot(aes(x=estimate, y=log(1/theta))) +
  stat_binhex(bins=100) +
  xlab("mu") +
  ylab("log(1/theta)") +
  scale_fill_viridis_c(option="C") +
  geom_smooth(method="loess", method.args=list(degree=1), se=FALSE) +
  theme_bw() +
  theme(aspect.ratio = 1)

#assign the loess trend to object
Loess.Trend <- loess(log(1/theta)~estimate, data = Prop_Field, degree = 1)

#Add a column for the residual, which we will call 'dispersion'
Prop_Field <- Prop_Field %>%
  ungroup() %>%
  drop_na() %>%
  mutate(Dispersion = residuals(Loess.Trend))

#Plot the residuals, with a horizontal blue line as a way to visually communicate that we are plotting the residuals from the loess trend
Prop_Field %>%
  ggplot(aes(x=estimate, y=Dispersion)) +
  stat_binhex(bins=100) +
  xlab("mu") +
  ylab("Dispersion") +
  geom_hline(yintercept = 0, color="blue") +
  scale_fill_viridis_c(option="C") +
  theme_bw() +
  theme(aspect.ratio = 1)

write.csv(Prop_Field, quote=F, row.names=F, file="Prop_Field_Disp.csv")


