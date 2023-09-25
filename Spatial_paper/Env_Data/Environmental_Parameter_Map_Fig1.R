#Script makes the subplots shown in Figures 1 and S1. Subplots were assembled into figure outside of R. 

rm(list = ls())
library(stringr)
library(ggOceanMapsData)
library(ggOceanMaps)
library(viridis)
library(dplyr)
#throughout script replace "yourpath" with path to file
data <- read.csv("yourpath/Zooplankton_Density_LMG1901.csv", header=TRUE)

CTD_files <- list.files(path = "yourpath", pattern = "*.asc")
rosetta = read.csv("yourpath/CTD_Rosetta_Stone.csv")
#get rosetta cruise tow into same format
names(rosetta)[1] = "CruiseTow"
rosetta[,1] = rosetta[,1] + 270000
data = merge(data, rosetta, by = "CruiseTow")

# Note: the function is finicky...need to have only the desired *.asc files in the directory. No other *.asc files or *.asc.gz files

read_CTD_fun = function(cast){
  ctd_ids = as.numeric( tools::file_path_sans_ext( str_split(CTD_files, pattern="_", n=2, simplify=TRUE)[,2] ) )
  
  CTD <- read.table(CTD_files[ctd_ids == cast], sep="", header=TRUE, row.names=NULL)
  names(CTD)[1] = "Depth"
  max_temp = max(CTD$T090C)
  mean_temp = mean(CTD$T090C)
  min_temp = min(CTD$T090C)
  max_salinity = max(CTD$Sal11)
  mean_salinity = mean(CTD$Sal11)
  min_salinity = min(CTD$Sal11)
  max_chl = max(CTD$FlECO.AFL)
  shallow = CTD[which(CTD$DepSM <200),]
  int_chl = sum(shallow$FlECO.AFL)
  upper_temp = mean(shallow$T090C)
  upper_max = max(shallow$T090C)
  max_depth = max(shallow$DepSM)
  
  return(c(max_chl, int_chl, upper_temp, upper_max, max_depth))
}

CTD_files <- list.files(path = "C:/Users/atarrant/Documents/Paper_Ant_Spatial/Env_Data", pattern = "*.asc")
setwd("C:/Users/atarrant/Documents/Paper_Ant_Spatial/Env_Data")
ctd_data = sapply(data$Cast, read_CTD_fun)

row.names(ctd_data) = c("max_chl", "int_chl", "upper_temp", "upper_max_temp", "shallow lim")

data = cbind(data, t(ctd_data))
names(data)[names(data) == 'CopepodaNum..num.1000m?.'] <- 'Copepods_per_1000_m3'
names(data)[names(data) == 'CopepodaVol..ml.1000m3.'] <- 'Copepods_ml_per_m3'

#output datafile
#write.csv(data, "Station_summary.csv")

data$Lat <- -1*(data$LatitudeStart..?.)
data$Long <- -1*(data$LongitudeStart..?.)
sampled <- filter(data, Sampled == 2)
not <- filter(data, Sampled == 1)
counts <- filter(data, C_pro != "NaN")

a<-basemap(limits = c(-2.8e6, -2e6, 4e5, 1.4e6), projection.grid = FALSE, shapefiles = "Antarctic", bathymetry=TRUE, bathy.style = "contour_blues", legends = FALSE, grid.col=NA)+
  annotation_scale()+
  labs(x=NULL, y=NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_spatial_point(data=data, aes(x=Long, y=Lat, fill=Copepods_ml_per_m3), shape=21, size = 5, color="black")+
  labs(fill=expression("Copepods, ml/1000" ~m^3))+
  theme(legend.position = c(0.83, 0.85))+
  scale_fill_viridis(option="D")

b<-basemap(limits = c(-2.8e6, -2e6, 4e5, 1.4e6), projection.grid = FALSE, shapefiles = "Antarctic", bathymetry=TRUE, bathy.style = "contour_blues", legends = FALSE, grid.col=NA)+
  annotation_scale()+
  labs(x=NULL, y=NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_spatial_point(data=counts, aes(x=Long, y=Lat, fill=C_acu), shape=21, size=5, color="black")+
  labs(fill=expression("Calanoides / 1000"~m^3))+
  theme(legend.position = c(0.83, 0.85))+
  scale_fill_viridis(option="D")

c<-basemap(limits = c(-2.8e6, -2e6, 4e5, 1.4e6), projection.grid = FALSE, shapefiles = "Antarctic", bathymetry=TRUE, bathy.style = "contour_blues", legends = FALSE, grid.col=NA)+
  annotation_scale()+
  labs(x=NULL, y=NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_spatial_point(data=counts, aes(x=Long, y=Lat, fill=C_pro), shape=21, size=5, color="black")+
  labs(fill=expression("Calanus / 1000" ~m^3))+
  theme(legend.position = c(0.83, 0.85))+
  scale_fill_viridis(option="D")

data$int_chl <- data$int_chl/200

d<- basemap(limits = c(-2.8e6, -2e6, 4e5, 1.4e6), projection.grid = FALSE, shapefiles = "Antarctic", bathymetry=TRUE, bathy.style = "contour_blues", legends = FALSE, grid.col=NA)+
  annotation_scale()+
  labs(x=NULL, y=NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_spatial_point(data=data, aes(x=Long, y=Lat, fill=int_chl), shape=21, size=5, color="black")+
  labs(fill=expression("Chlorophyll mg/" ~m^3))+
  theme(legend.position = c(0.83, 0.85))+
  scale_fill_viridis(option="D")

#Temperature figures for supplement
e<- basemap(limits = c(-2.8e6, -2e6, 4e5, 1.4e6), projection.grid = FALSE, shapefiles = "Antarctic", bathymetry=TRUE, bathy.style = "contour_blues", legends = FALSE, grid.col=NA)+
  annotation_scale()+
  labs(x=NULL, y=NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_spatial_point(data=data, aes(x=Long, y=Lat, fill=upper_max_temp), shape=21, size=5, color="black")+
  labs(fill=expression("Max Upper Temp, C"))+
  theme(legend.position = c(0.83, 0.85))+
  scale_fill_viridis(option="D")

f<- basemap(limits = c(-2.8e6, -2e6, 4e5, 1.4e6), projection.grid = FALSE, shapefiles = "Antarctic", bathymetry=TRUE, bathy.style = "contour_blues", legends = FALSE, grid.col=NA)+
  annotation_scale()+
  labs(x=NULL, y=NULL)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_spatial_point(data=data, aes(x=Long, y=Lat, fill=upper_temp), shape=21, size=5, color="black")+
  labs(fill=expression("Avg. Temp 0-200 m, C"))+
  theme(legend.position = c(0.83, 0.85))+
  scale_fill_viridis(option="D")

