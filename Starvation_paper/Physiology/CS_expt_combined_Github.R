# Read in data. 
rm(list = ls())
#in () add the directory where the data are located
setwd()
Cp_Expt <- read.csv("CS_Cp_Experiment.csv", header=TRUE)
Cp_Expt$Time = factor(Cp_Expt$Time)
Cp_Expt$Group = factor(Cp_Expt$Group)
Ca_Expt <- read.csv("CS_Ca_Experiment.csv", header=TRUE)
Ca_Expt$Time = factor(Ca_Expt$Time)
Ca_Expt$Group = factor(Ca_Expt$Group)
# Adjust values for 0.3 ml homogenization volume
Ca_Expt$CS_weight <- (Ca_Expt$CS/Ca_Expt$Wet_Weight)*.3
Cp_Expt$CS_weight <- (Cp_Expt$CorrCS/Cp_Expt$Wet_Weight)*.3

# "userfriendlyscience" has been replaced by "ufs"
library(ufs)
library(car)
library(ggplot2)
library(dplyr)
library(gridExtra)

#Propinquus data
#Remove time 0 from Feeding experiment to have two-way ANOVA design
Cp_Expt_2T <- Cp_Expt[ which(Cp_Expt$Time!="0"),]


#Normalized to weight
m <- aov(CS_weight ~ Time * Group, data = Cp_Expt_2T)
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
leveneTest(CS_weight ~ Time * Group, data = Cp_Expt_2T, center=mean) # Levene's test, violates assumption
# Can adjust for heteroscedasticity using ANOVA within 'car'
mod <- lm(CS_weight ~ Time * Group, data=Cp_Expt_2T, contrasts=list(Time=contr.sum, Group=contr.sum))
Anova(mod, type=3, white.adjust=T)  # note use of contr.sum in call to lm()
# Significant effect of Group only. Not that we get the same (robust) result without the white adjustment.


#Acutus data
Ca_Expt_2T <- Ca_Expt[ which(Ca_Expt$Time!="0" & Ca_Expt$Group!="Refed"),]

#Normalized to weight
m <- aov(CS_weight ~ Time * Group, data = Ca_Expt_2T)
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
leveneTest(CS_weight ~ Time * Group, data = Ca_Expt_2T, center=mean) # passes Levene's test
summary(m)
# Significant effect of time and group.


Cp_plot <-Cp_Expt_2T %>%
  ggplot(aes(x=Time, y=CS_weight, color=Group)) +
  geom_boxplot() +
  ylab("Citrate Synthase") +
  xlab("Time") +
  #coord_fixed() +
  theme_minimal()

Ca_plot <- Ca_Expt_2T %>%
  ggplot(aes(x=Time, y=CS_weight, color=Group)) +
  geom_boxplot() +
  ylab("Citrate Synthase") +
  xlab("Time") +
  #coord_fixed() +
  theme_minimal()


both <- grid.arrange(Cp_plot, Ca_plot, ncol = 2, nrow=2)

ggsave(plot = both, file= "combo_plot_CS_weight.pdf")
ggsave(plot = both, file = "combo_plot_CS_weight.png")
