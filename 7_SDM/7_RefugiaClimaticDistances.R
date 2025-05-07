
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

sample.info = read.table("../0_SampleInfo.txt", header = T)
info = read.csv("../../InfOwlsEurope_June2020.csv", header = T)
subset.info = info[match(sample.info$ID, info$LibPapers), ]
ITPT.GPS = data.frame(long = subset.info$Longitude_decimal[sample.info$population == "Portugal" | sample.info$population == "Italy"], 
                      lat = subset.info$Latitude_decimal[sample.info$population == "Portugal" | sample.info$population == "Italy"])
ITPT.ID = substr(c(sample.info$ID[grepl("IT", sample.info$ID)], sample.info$ID[grepl("PT", sample.info$ID)]), 1, 2)

predictors.past.ITPT = raster::extract(PastPredictors, ITPT.GPS)
predictors.present.ITPT = raster::extract(CurrentPredictors, ITPT.GPS)

PCA.ITPT.past = predict(PCA, predictors.past.ITPT)
PCA.ITPT.present = predict(PCA, predictors.present.ITPT)
rownames(PCA.ITPT.present) = rownames(PCA.ITPT.past) = ITPT.ID

###############################

ITPT.ID = paste0(ITPT.ID, "0", rep(1:9, 2))

pairwise.eucD.past = matrix(NA, nrow(PCA.ITPT.past), ncol = nrow(PCA.ITPT.past))
colnames(pairwise.eucD.past) = rownames(pairwise.eucD.past) = ITPT.ID

for (i in 1:nrow(PCA.ITPT.past)) {
  row.i = PCA.ITPT.past[i, 1:3]
  for (j in 1:nrow(PCA.ITPT.past)) {
    row.j = PCA.ITPT.past[j, 1:3]
    eucD = as.numeric(sqrt(((row.i[1] - row.j[1])^2) +
                             ((row.i[2] - row.j[2])^2) +
                             ((row.i[3] - row.j[3])^2) ))
    pairwise.eucD.past[i, j] = eucD
  }
}

pairwise.eucD.past = pairwise.eucD.past[grepl("PT", rownames(pairwise.eucD.past)), 
                                        grepl("IT", colnames(pairwise.eucD.past))]
pairwise.eucD.past = pairwise.eucD.past[lower.tri(pairwise.eucD.past, diag = F)] 

pairwise.eucD.present = matrix(NA, nrow = nrow(PCA.ITPT.present), ncol = nrow(PCA.ITPT.present))
colnames(pairwise.eucD.present) = rownames(pairwise.eucD.present) = ITPT.ID

for (i in 1:nrow(PCA.ITPT.present)) {
  row.i = PCA.ITPT.present[i, 1:3]
  for (j in 1:nrow(PCA.ITPT.present)) {
    row.j = PCA.ITPT.present[j, 1:3]
    eucD = as.numeric(sqrt(((row.i[1] - row.j[1])^2) +
                             ((row.i[2] - row.j[2])^2) +
                             ((row.i[3] - row.j[3])^2) ))
    pairwise.eucD.present[i, j] = eucD
  }
}

pairwise.eucD.present = pairwise.eucD.present[grepl("PT", rownames(pairwise.eucD.present)), 
                                              grepl("IT", colnames(pairwise.eucD.present))]
pairwise.eucD.present = pairwise.eucD.present[lower.tri(pairwise.eucD.present, diag = F)] 
