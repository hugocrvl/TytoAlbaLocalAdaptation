
###############################

########### LOCALLY ###########

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../../0_Functions.R", chdir = T)

###############################

library("rbioclim")
library("rgl")

###############################

sample.info = read.table("", header = T) # Individual and population IDs
info = read.csv("", header = T) # Dataset metadata
subset.info = info[match(sample.info$ID, info$LibPapers), ]

CurrentEnv = rbioclim::getData("worldclim", var = "bio", res = 5, 
                               path = "")
bioclim.owl.sampling = data.frame(raster::extract(CurrentEnv,
                                                  data.frame(x = subset.info$Longitude_decimal,
                                                             y = subset.info$Latitude_decimal)))
RDA.explanatory = bioclim.owl.sampling[, c("bio2","bio6","bio7","bio8", "bio15","bio17", "bio19")] # See 7_SDM scripts for variable selection
write.table(RDA.explanatory, "", quote = F, row.names = F, col.names = T)

