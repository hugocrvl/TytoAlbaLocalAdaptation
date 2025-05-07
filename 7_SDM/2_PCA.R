
###############################

########### LOCALLY ###########

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../0_Functions.R", chdir = T)

###############################

# Loading packages
library("dismo")
library("rgeos")
# library("rgdal")
library("mapplots")
library("maps")
library("mapdata")
library("maptools")
library("jsonlite")
# library("rJava")
library("ENMeval")
library("rbioclim")
library("splancs")
library("rgl")
library("alphashape3d")
library("factoextra")

###############################

# Area around Europe and North-African regions
area = as(raster::extent(-20, 46, 22, 68), "SpatialPolygons")

# Load bioclimatic variables from the past and now
PastEnv = rbioclim::getData("worldclim_past", var = "bio", res = 5,
                            past = "lgm", path = "0_WorldClim")
CurrentEnv = rbioclim::getData("worldclim", var = "bio", res = 5, 
                               path = "0_WorldClim")
names(PastEnv) = names(CurrentEnv)

# Crop around Europe And North-African regions
CurrentEnv.crop = crop(CurrentEnv, area)
PastEnv.crop = crop(PastEnv, area)

CurrentPredictors = subset(CurrentEnv.crop, c("bio2","bio6","bio7","bio8", "bio15","bio17", "bio19"))
PastPredictors = subset(PastEnv.crop, c("bio2","bio6","bio7","bio8", "bio15","bio17", "bio19"))

################################################################################

mean.predict.past = raster("0_Rasters/Mean_Predict_Past_MESS.grd")
mean.predict.present = raster("0_Rasters/Mean_Predict_Now_MESS.grd")

#####                                                            #####
#### Apply a PCA on the every continental pixel of the study area ####
#####                                                            #####

# Polygon around all continental pixel from the past and present
continent.past = mean.predict.past
continent.past[continent.past == 0] = 1
poly.continent.past = rasterToPolygons(continent.past, dissolve = T)

continent.now = mean.predict.present
continent.now[continent.now == 0] = 1
poly.continent.present = rasterToPolygons(continent.now, dissolve = T)

# Extraction of the climatic variables inside these two polygons and merging of all values
bioclim.continent.past = data.frame(raster::extract(PastPredictors, poly.continent.past))
bioclim.continent.present = data.frame(raster::extract(CurrentPredictors, poly.continent.present))
overall.bioclim.continent = rbind(bioclim.continent.past, bioclim.continent.present)

# PCA on all these climatic variables, and creation of a vector with explained variance per axis
PCA = princomp(overall.bioclim.continent)

summary(PCA)

load = loadings(PCA)

principal.components = data.frame(PCA$scores)
index.past.present = c(rep("past", nrow(bioclim.continent.past)), rep("present", nrow(bioclim.continent.present)))

# Distinction between the unsuitable and suitable pixels from continental area from the past and present
# Useful when I'll select the suitable pixel only inside the PCA space
binary.suitable.past = data.frame(raster::extract(mean.predict.past, poly.continent.past))
binary.suitable.present = data.frame(raster::extract(mean.predict.present, poly.continent.now))
colnames(binary.suitable.past) = colnames(binary.suitable.present) = "index"

binary.suitable.both = rbind(binary.suitable.past, binary.suitable.present)
output.PCA = cbind(principal.components, index.past.present, binary.suitable = binary.suitable.both$index)

write.table(output.PCA, "2.1_PCA.txt", quote = F, col.names = T, row.names = F)
