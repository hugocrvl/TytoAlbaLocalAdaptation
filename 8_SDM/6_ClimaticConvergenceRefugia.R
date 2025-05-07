
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
library("scales")

###############################

area = as(raster::extent(-20, 46, 22, 68), "SpatialPolygons")

CurrentEnv = rbioclim::getData("worldclim", var = "bio", res = 5, 
                               path = "0_WorldClim")
PastEnv = rbioclim::getData("worldclim_past", var = "bio", res = 5,
                            past = "lgm", path = "0_WorldClim")
names(PastEnv) = names(CurrentEnv)

CurrentEnv.crop = crop(CurrentEnv, area)
PastEnv.crop = crop(PastEnv, area)

CurrentPredictors = subset(CurrentEnv.crop, c("bio2","bio6","bio7","bio8", "bio15","bio17", "bio19"))
PastPredictors = subset(PastEnv.crop, c("bio2","bio6","bio7","bio8", "bio15","bio17", "bio19"))

mean.predict.past = raster("0_Rasters/Mean_Predict_Past_MESS.grd")
mean.predict.present = raster("0_Rasters/Mean_Predict_Now_MESS.grd")

# The PCA needs to be done again to be able to re project the bioclimatic
# variables from the common suitable area
continent.past = mean.predict.past
continent.past[continent.past == 0] = 1
poly.continent.past = rasterToPolygons(continent.past, dissolve = T)

continent.present = mean.predict.present
continent.present[continent.present != 0] = NA
polygon.continent.present = rasterToPolygons(continent.present, dissolve = T)

bioclim.continent.past = data.frame(raster::extract(PastPredictors, poly.continent.past))
bioclim.continent.present = data.frame(raster::extract(CurrentPredictors, polygon.continent.present))
overall.bioclim.continent = rbind(bioclim.continent.past, bioclim.continent.present)
PCA = princomp(overall.bioclim.continent)

# Extraction of the bioclimatic variables from the common suitable area
mean.predict.common = mean.predict.past + mean.predict.present
mean.predict.common[mean.predict.common != 2] = NA
polygon.suitable.common = rasterToPolygons(mean.predict.common, dissolve = T)

bioclim.common.past = data.frame(raster::extract(PastPredictors, polygon.suitable.common))
bioclim.common.present = data.frame(raster::extract(CurrentPredictors, polygon.suitable.common))

# Projection of theses variables in the PCA space
projection.past.common = data.frame((predict(PCA, bioclim.common.past)))
projection.present.common = data.frame((predict(PCA, bioclim.common.present)))

distance.past.present = c()
for (i in 1:nrow(projection.past.common)) {
  
  row.past = projection.past.common[i, 1:3]
  row.present = projection.present.common[i, 1:3]
  distance = as.numeric(sqrt(((row.present[1] - row.past[1])^2)) + 
                          (((row.present[2] - row.past[2])^2)) + 
                          (((row.present[3] - row.past[3])^2)))
  distance.past.present[i] = distance
  
}

id.cell.common = extract(mean.predict.present, polygon.suitable.common, cellnumbers = TRUE)[[1]]
row.column.common = raster::rowColFromCell(mean.predict.present, id.cell.common[, 1])

# distance.past.present.rescaled = distance.past.present / max(distance.past.present)
n = length(distance.past.present)
distance.past.present.rescaled = rescale(distance.past.present, c(0,1))
palette = colorRampPalette(c("#046235", "#046235", "#25F892", "#FFFF66"))(n)
colour.distance.common = palette[as.factor(distance.past.present.rescaled)]

plot(distance.past.present, col = colour.distance.common, pch = 20, ylab = "Distance")

mean.predict.common.unsuitable = mean.predict.past + mean.predict.present
mean.predict.common.unsuitable[mean.predict.common.unsuitable == 2] = NA
mean.predict.common.unsuitable[mean.predict.common.unsuitable == 1] = 0
polygon.unsuitable.common = rasterToPolygons(mean.predict.common.unsuitable, dissolve = T)

colour.distance.common.factor = as.factor(colour.distance.common)
r = mean.predict.present
r[] = NA
r[row.column.common] = as.numeric(colour.distance.common.factor)

pdf("6.1_EuropeChange.pdf", height = 15, width = 15)

plot(polygon.unsuitable.common, density = 15, angle = 0, col = "gray90")
plot(polygon.unsuitable.common, density = 15, angle = 90, col = "gray90", add = T)

plot(r, col = levels(colour.distance.common.factor), add = TRUE, legend = F)
plot(coast.area, add = TRUE, col = add.alpha("#EFEFEF", alpha = 0.4),
     legend = F, axes = F)

dev.off()