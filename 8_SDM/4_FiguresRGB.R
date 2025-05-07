
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

library("rworldmap")
library("dplyr")
library("ggplot2")
library("raster")
library("sf")

library("extrafont")
loadfonts()

# library("RCGLS")
# library("geosphere")
# library(gpclib)
# library("GSIF")
# library("openeo")
# library("tibble")

###############################

# Coordinates from our samples

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

populations = substr(sample.info$ID, 1, 2)
list.colour.bg = colour.population[as.factor(populations)]
list.pch = c(21:25, 21:24)[as.factor(populations)]

###############################

principal.components = read.table("2.1_PCA.txt", header = T)
inside.3d = read.table("3.1_InsidePolygon.txt", header = F)$V1

principal.components.present = principal.components[principal.components$index.past.present == "present", ]
principal.components.present.suitable = principal.components.present[principal.components.present$binary.suitable == 1, ]

colour.inside.outside.3d = choose.colours(principal.components.present.suitable$Comp.1,
                                                 principal.components.present.suitable$Comp.2,
                                                 principal.components.present.suitable$Comp.3)
colour.inside.outside.3d[inside.3d] = "gray40"
colour.inside.outside.3d.factor = as.factor(colour.inside.outside.3d)

#
plot3d(principal.components.present.suitable, col = colour.inside.outside.3d)
#

mean.predict.present = raster("0_Rasters/Mean_Predict_Now_MESS.grd")

continent.present = mean.predict.present
continent.present[continent.present != 0] = NA
polygon.continent.present = rasterToPolygons(continent.present, dissolve = T)

mean.predict.present.only.suitable = mean.predict.present
mean.predict.present.only.suitable[mean.predict.present.only.suitable != 1] = NA
polygon.suitable.present = rasterToPolygons(mean.predict.present.only.suitable, dissolve = T)
id.cell.present = extract(mean.predict.present, polygon.suitable.present, cellnumbers = TRUE)[[1]]
row.column.present = raster::rowColFromCell(mean.predict.present, id.cell.present[, 1])
mean.predict.present.only.suitable[row.column.present] = as.numeric(colour.inside.outside.3d.factor)

binary.inside.outside = rep(0, nrow(row.column.present))
binary.inside.outside[inside.3d] = 1
temp.inside.polygon = mean.predict.present
temp.inside.polygon[temp.inside.polygon != 1] = NA
temp.inside.polygon[row.column.present] = binary.inside.outside
temp.inside.polygon[temp.inside.polygon != 1] = NA
inside.polygon.3D = rasterToPolygons(temp.inside.polygon, dissolve = TRUE)

png("8.1_EuropeRGB.png", height = 4500, width = 4500, res = 300)
# pdf("8.1_EuropeRGB.pdf", height = 15, width = 15)

plot(polygon.continent.present, density = 15, angle = 0, col = "gray90")
plot(polygon.continent.present, density = 15, angle = 90, col = "gray90", add = T)

plot(mean.predict.present.only.suitable, col = levels(colour.inside.outside.3d.factor),
     add = TRUE, legend = F)
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


png("8.4_LegendPopulations.png", height = 1500, width = 4500, res = 300)
# pdf("8.4_LegendPopulations.pdf", height = 5, width = 15)

plot(1:10, 1:10, type = 'n', axes = F, xlab = '', ylab = '', family = "Papyrus")
legend("center", legend = unique(sample.info$population), col = colour.population, 
       pt.bg = colour.population, pch = c(21:25, 21:24), 
       ncol = 5, cex = 2, pt.cex = 2)
legend("center", legend = unique(sample.info$population), col = c(rep(NA, 5), rep("#000000", 4)),
       pch = 20, text.col = add.alpha("#000000", alpha = 0), bg = NA, 
       ncol = 5, cex = 2, pt.cex = 1)

dev.off()

png("8.4_LegendPopulationsDRAFT.png", height = 3000, width = 1500, res = 300)

plot(1:10, 1:10, type = 'n', axes = F, xlab = '', ylab = '', family = "Papyrus")
legend("center", legend = unique(sample.info$population), col = colour.population, 
       pt.bg = colour.population, pch = c(21:25, 21:24), 
       ncol = 1, cex = 2, pt.cex = 2)
legend("center", legend = unique(sample.info$population), col = c(rep(NA, 5), rep("#000000", 4)),
       pch = 20, text.col = add.alpha("#000000", alpha = 0), bg = NA, 
       ncol = 1, cex = 2, pt.cex = 1)

dev.off()

png("../../../PhD/Conferences/Biology25/LegendPopulations.png", height = 3000, width = 1500, res = 300,
    bg = "transparent")

plot(1:10, 1:10, type = 'n', axes = F, xlab = '', ylab = '', family = "Papyrus")
legend("center", legend = unique(sample.info$population), col = colour.population, 
       pt.bg = colour.population, pch = c(21:25, 21:24), 
       ncol = 1, cex = 2, pt.cex = 2, box.col = NA)
legend("center", legend = unique(sample.info$population), col = c(rep(NA, 5), rep("#000000", 4)),
       pch = 20, text.col = add.alpha("#000000", alpha = 0), bg = NA, 
       ncol = 1, cex = 2, pt.cex = 1, box.col = NA)

dev.off()

###############################

worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id

world.df <- world.points[,c("long","lat","group", "region")]

test = data.frame(c(40,-21,-21,40),
                  c(25,25,50,50),
                  c("Zoom","Zoom","Zoom","Zoom"),
                  c("Zoom","Zoom","Zoom","Zoom"))
colnames(test) = c("long", "lat", "group", "region")

test = data.frame(c(46,-20,-20,46),
                  c(22,22,68,68),
                  c("Zoom","Zoom","Zoom","Zoom"),
                  c("Zoom","Zoom","Zoom","Zoom"))
colnames(test) = c("long", "lat", "group", "region")

worldmap = ggplot() +
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) +
  geom_polygon(data = test, aes(x = long, y = lat, group = region, col = "red",  alpha = 0.1)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  coord_map("ortho", orientation=c(25, 10, 0))
worldmap

worldmap

pdf("8.5_MapWorld.pdf", height = 15, width = 15)
png("8.5_MapWorld.png", height = 4500, width = 4500, res = 300)

worldmap

dev.off()
