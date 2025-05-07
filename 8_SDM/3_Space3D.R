
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

principal.components = read.table("2.1_PCA.txt", header = T)
explained.variance.percentage = read.table("2.2_ExplainedVariance.txt", header = F)$V1
principal.components.past = principal.components[principal.components$index.past.present == "past", ]
principal.components.present = principal.components[principal.components$index.past.present == "present", ]

###############################

principal.components.normalised = matrix(ncol = 3, 
                                         nrow = nrow(principal.components))
for (i in 1:3) {
  
  principal.components.normalised[, i] = (principal.components[, i] - 
                                            mean(principal.components[, i])) / sd(principal.components[, i])
  
}

principal.components.normalised.past = as.matrix(cbind(principal.components.normalised[principal.components$index.past.present == "past", 1],
                                                       principal.components.normalised[principal.components$index.past.present == "past", 2],
                                                       principal.components.normalised[principal.components$index.past.present == "past", 3]))
principal.components.normalised.past.suitable = principal.components.normalised.past[principal.components.past$binary.suitable == 1, ]


principal.components.normalised.present = as.matrix(cbind(principal.components.normalised[principal.components$index.past.present == "present", 1],
                                                          principal.components.normalised[principal.components$index.past.present == "present", 2],
                                                          principal.components.normalised[principal.components$index.past.present == "present", 3]))
principal.components.normalised.present.suitable = principal.components.normalised.present[principal.components.present$binary.suitable == 1, ]

ashape.obj.past = ashape3d(as.matrix(principal.components.normalised.past.suitable), alpha = 0.5)
inside.3d = inashape3d(ashape.obj.past, points = principal.components.normalised.present.suitable)
table(inside.3d)
# FALSE  TRUE 
# 56044 19224

write.table(inside.3d, "3.1_InsidePolygon.txt", quote = F, row.names = F,
            col.names = F)
