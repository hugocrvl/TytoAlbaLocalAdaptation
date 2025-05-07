
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
library("scales")

###############################

set.seed(0202)
sample.info = read.table("../0_SampleInfo.txt", header = T)
info = read.csv("../../InfOwlsEurope_June2020.csv", header = T)
subset.info = info[match(sample.info$ID, info$LibPapers), ]
coordinates = data.frame(lon = subset.info$Longitude_decimal, 
                         lat = subset.info$Latitude_decimal)
lon = jitter(coordinates$lon, amount = 0.425)
lat = jitter(coordinates$lat, amount = 0.42)

area = as(raster::extent(-20, 46, 22, 68), "SpatialPolygons")
coast = shapefile("0_Coast/ne_10m_coastline.shp")
proj4string(coast) = crs(area)
coast.area = gIntersection(coast, area)

mean.predict.present = raster("0_Rasters/Mean_Predict_Now_MESS.grd")
inside.3d = read.table("3.1_InsidePolygon.txt", header = F)$V1

continent.present = mean.predict.present
continent.present[continent.present != 0] = NA
polygon.continent.present = rasterToPolygons(continent.present, dissolve = T)

principal.components = read.table("2.1_PCA.txt", header = T)
principal.components.present = principal.components[principal.components$index.past.present == "present", ]
principal.components.present.suitable = principal.components.present[principal.components.present$binary.suitable == 1, ]

mean.predict.present.only.suitable = mean.predict.present
mean.predict.present.only.suitable[mean.predict.present.only.suitable != 1] = NA
polygon.suitable.present = rasterToPolygons(mean.predict.present.only.suitable, dissolve = T)

mean.predict.present.only.suitable.PC1 = mean.predict.present.only.suitable.PC2 = 
  mean.predict.present.only.suitable.PC3 = mean.predict.present.only.suitable

id.cell.present = extract(mean.predict.present, polygon.suitable.present, cellnumbers = TRUE)[[1]] # Extraction of the present suitable pixels
row.column.present = raster::rowColFromCell(mean.predict.present, id.cell.present[, 1]) # Get the pixel index

binary.inside.outside = rep(0, nrow(row.column.present))
binary.inside.outside[inside.3d] = 1
temp.inside.polygon = mean.predict.present
temp.inside.polygon[temp.inside.polygon != 1] = NA
temp.inside.polygon[row.column.present] = binary.inside.outside
temp.inside.polygon[temp.inside.polygon != 1] = NA
inside.polygon.3D = rasterToPolygons(temp.inside.polygon, dissolve = TRUE)

n = nrow(principal.components.present.suitable)
palette = colorRampPalette(c("#000000", "#FFFFFF"))(n)

rescale.PC1 = rescale(principal.components.present.suitable$Comp.1, c(0,1))
colour.PC1 = palette[as.factor(rescale.PC1)]
colour.PC1.factor = as.factor(colour.PC1)
mean.predict.present.only.suitable.PC1[row.column.present] = as.numeric(colour.PC1.factor)

png("9.1_ContributionPC1Europe.png", height = 4500, width = 4500, res = 300)
# pdf("9.1_ContributionPC1Europe.pdf", height = 15, width = 15)

plot(polygon.continent.present, density = 15, angle = 0, col = "gray90")
plot(polygon.continent.present, density = 15, angle = 90, col = "gray90", add = T)

plot(mean.predict.present.only.suitable.PC1, col = levels(colour.PC1.factor),
     add = T, legend = F, axes = F)
plot(coast.area, add = TRUE, col = add.alpha("#EFEFEF", alpha = 0.4),
     legend = F, axes = F)
plot(inside.polygon.3D, border = "#EFEFEF", add = TRUE)

for (i in 1:length(populations)) {
  
  points(lon[i], lat[i], col = list.colour.bg[i], bg = list.colour.bg[i], pch = list.pch[i],
         cex = 2.25)
  
  if (populations[i] %in% unique(populations)[6:9]) {
    
    points(lon[i], lat[i], pch = 20, col = "#000000", cex = 1)
    
  }
  
}

dev.off()

rescale.PC2 = rescale(principal.components.present.suitable$Comp.2, c(0,1))
colour.PC2 = palette[as.factor(rescale.PC2)]
colour.PC2.factor = as.factor(colour.PC2)
mean.predict.present.only.suitable.PC2[row.column.present] = as.numeric(colour.PC2.factor)

png("9.2_ContributionPC2Europe.png", height = 4500, width = 4500, res = 300)
# pdf("9.2_ContributionPC2Europe.pdf", height = 15, width = 15)

plot(polygon.continent.present, density = 15, angle = 0, col = "gray90")
plot(polygon.continent.present, density = 15, angle = 90, col = "gray90", add = T)

plot(mean.predict.present.only.suitable.PC2, col = levels(colour.PC2.factor),
     add = T, legend = F, axes = F)
plot(coast.area, add = TRUE, col = add.alpha("#EFEFEF", alpha = 0.4),
     legend = F, axes = F)
plot(inside.polygon.3D, border = "#EFEFEF", add = TRUE)

for (i in 1:length(populations)) {
  
  points(lon[i], lat[i], col = list.colour.bg[i], bg = list.colour.bg[i], pch = list.pch[i],
         cex = 2.25)
  
  if (populations[i] %in% unique(populations)[6:9]) {
    
    points(lon[i], lat[i], pch = 20, col = "#000000", cex = 1)
    
  }
  
}

dev.off()

rescale.PC3 = rescale(principal.components.present.suitable$Comp.3, c(0,1))
colour.PC3 = palette[as.factor(rescale.PC3)]
colour.PC3.factor = as.factor(colour.PC3)
mean.predict.present.only.suitable.PC3[row.column.present] = as.numeric(colour.PC3.factor)

png("9.3_ContributionPC3Europe.png", height = 4500, width = 4500, res = 300)
# pdf("9.3_ContributionPC3Europe.pdf", height = 15, width = 15)

plot(polygon.continent.present, density = 15, angle = 0, col = "gray90")
plot(polygon.continent.present, density = 15, angle = 90, col = "gray90", add = T)

plot(mean.predict.present.only.suitable.PC3, col = levels(colour.PC3.factor),
     add = T, legend = F, axes = F)
plot(coast.area, add = TRUE, col = add.alpha("#EFEFEF", alpha = 0.4),
     legend = F, axes = F)
plot(inside.polygon.3D, border = "#EFEFEF", add = TRUE)

for (i in 1:length(populations)) {
  
  points(lon[i], lat[i], col = list.colour.bg[i], bg = list.colour.bg[i], pch = list.pch[i],
         cex = 2.25)
  
  if (populations[i] %in% unique(populations)[6:9]) {
    
    points(lon[i], lat[i], pch = 20, col = "#000000", cex = 1)
    
  }
  
}

dev.off()
