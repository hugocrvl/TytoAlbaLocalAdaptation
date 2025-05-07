
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
library("rJava")
library("ENMeval")
library("rbioclim")
library("splancs")
library("rgl")
library("alphashape3d")
library("ade4")

###############################

# Load barn owl distribution area
BarnOwlRange = shapefile("../../environmental_niche/RunSDM_Hugo/IUCN/data_0.shp")
# Crop around Europe and North-African regions
area = as(raster::extent(-20, 46, 22, 68), "SpatialPolygons")
proj4string(area) = crs(BarnOwlRange)
RangeArea = gIntersection(BarnOwlRange, area)

# Load bioclimatic variables from the past and now
CurrentEnv = rbioclim::getData("worldclim", var = "bio", res = 5, 
                               path = "0_WorldClim")
PastEnv = rbioclim::getData("worldclim_past", var = "bio", res = 5,
                            past = "lgm", path = "0_WorldClim")
names(PastEnv) = names(CurrentEnv)

# Crop around Europe And North-African regions
CurrentEnv.crop = crop(CurrentEnv, area)
PastEnv.crop = crop(PastEnv, area)

CurrentPredictors = subset(CurrentEnv.crop, c("bio2","bio6","bio7","bio8", "bio15","bio17", "bio19"))
PastPredictors = subset(PastEnv.crop, c("bio2","bio6","bio7","bio8", "bio15","bio17", "bio19"))

# Subset the variables based on a correlation below 0.8
png("1.7_SelectClimaticVariables.png", width = 3000, height = 3000, res = 300)

pairs(CurrentPredictors)

dev.off()
# CurrentEnv.cor = data.frame(round(cor(CurrentEnv), 3))

coast = shapefile("0_Coast/ne_10m_coastline.shp")
proj4string(coast) = crs(area)
coast.area = gIntersection(coast, area)

################################################################################

#####                  #####
#### Run 1 MaxEnt model ####
#####                  #####

# Get sampling points
set.seed(0202)
Points = spsample(RangeArea, n=1000, "random")

# Choose the best parameters
Summary_Models = matrix(nrow = 3, ncol = 5)
Summary_Models = as.data.frame(Summary_Models)
Summary_Models_AUC = matrix(nrow = 3, ncol = 5)
Summary_Models_AUC = as.data.frame(Summary_Models_AUC)

rownames(Summary_Models) = c("-l", "-q","-h") #
colnames(Summary_Models) = as.character(c(1,2,3,4,5)) # Beta multiplier
rownames(Summary_Models_AUC) = c("-l", "-q","-h") #
colnames(Summary_Models_AUC) = as.character(c(1,2,3,4,5)) # Beta multiplier

get_nzc = function(df_lambda) {
  ncoeff = 0
  for (i in 1:length(df_lambda)) {
    x = as.numeric(strsplit(df_lambda[i], ",")[[1]][2])
    if (x != 0) {ncoeff = ncoeff + 1}
  }
  return(as.numeric(ncoeff))
}

for (Feature in c("-l", "-q","-h")){
  print(Feature)
  for (Beta in c(1,2,3,4,5)){
    print(Beta)
    myModel = maxent(x=CurrentPredictors, p=Points, args=c(Feature, paste0("betamultiplier=",Beta)))
    # Default : max 500 iterations, 1 replicate (1 model therefore one partition)
    Prediction = predict(myModel, CurrentPredictors)
    AIC = aic.maxent(p.occs = data.frame(Points), ncoefs = get_nzc(myModel@lambdas), p = Prediction)
    Summary_Models[Feature, as.character(Beta)]=AIC$AICc[1]
    Summary_Models_AUC[Feature, as.character(Beta)]=myModel@results["Training.AUC",]
  }
}
# 
# Summary_Models
#       1        2        3        4        5
# -l 13987.71 13981.23 13994.70 14032.05 14072.61
# -q 13943.59 13970.03 13994.54 14041.18 14069.88
# -h 14016.46 14082.10 14138.61 14185.25 14214.70
# MIN -> q 1, recommended by Warren & Seifert 2011

# Summary_Models_AUC
#       1      2      3      4      5
# -l 0.7909 0.7761 0.7745 0.7737 0.7754
# -q 0.7902 0.7812 0.7764 0.7762 0.7781
# -h 0.7635 0.7604 0.7565 0.7575 0.7516

# Choose the best model from different replicates
myMultipleModel = maxent(x = CurrentPredictors, p = Points,
                         args=c("randomtestpoints=25","randomseed=true",
                                "maximumiterations=5000",
                                "replicates=5", "replicatetype=subsample", 
                                "-q", "betamultiplier=1"))

myBestModel = myMultipleModel@models[[which(myMultipleModel@results["Training.AUC",]==max(myMultipleModel@results["Training.AUC",]))]]
# plot(myBestModel)

# Project the suitability now and during the LGM 
Prediction_now = predict(myBestModel, CurrentPredictors)
Prediction_PastEnv = predict(myBestModel, PastPredictors)

# Transformation into binary output
MTSS = myBestModel@results["Maximum.training.sensitivity.plus.specificity.Cloglog.threshold",] # At home
# MTSS = myBestModel@results["Maximum.training.sensitivity.plus.specificity.logistic.threshold", ] # Laptop

# Past
BinaryPrediction_Past = Prediction_PastEnv
BinaryPrediction_Past[BinaryPrediction_Past < MTSS] = 0
BinaryPrediction_Past[BinaryPrediction_Past >= MTSS] = 1

# Now
BinaryPrediction_Now = Prediction_now
BinaryPrediction_Now[BinaryPrediction_Now < MTSS] = 0
BinaryPrediction_Now[BinaryPrediction_Now >= MTSS] = 1

plot(BinaryPrediction_Past) ; plot(BinaryPrediction_Now)

################################################################################

#####                  #####
#### MESS investigation ####
#####                  #####

List_MESS_Past = list()
List_MESS_Curent = list()

i = 1
for (var in c("bio2","bio6","bio7","bio8", "bio15","bio17", "bio19")) {
  MESS = mess(PastEnv[[var]], extract(CurrentEnv[[var]], Points))
  List_MESS_Past[[i]] = MESS
  MESS = mess(CurrentEnv[[var]],extract(CurrentEnv[[var]], Points))
  List_MESS_Curent[[i]] = MESS
  i = i + 1
}

Stack_MESS_Past = raster::stack(List_MESS_Past)
Stack_MESS_Past_min = raster::stackApply(Stack_MESS_Past, 
                                         indices = length(c("bio2","bio6",
                                                            "bio7","bio8",
                                                            "bio15","bio17",
                                                            "bio19")), 
                                         fun = min)
values(Stack_MESS_Past_min)[values(Stack_MESS_Past_min) < 0] = 0
values(Stack_MESS_Past_min)[values(Stack_MESS_Past_min) > 0] = 1

################################################################################

#####                     #####
#### 100 replicated models ####
#####                     #####

List_Model = list()
List_Predict_Now = list()
List_Predict_Past = list()
MTSS_MyModel = c()
TestAUC = c()

List_Predict_Now_AUC = list()
List_Predict_Past_AUC = list()

set.seed(0202)
nmodel = 100
for(i in 1:nmodel){
  print(i)
  Points = spsample(RangeArea, n=1000, "random")
  MyModel =  maxent(x=CurrentPredictors, p=Points, args=c("randomtestpoints=25","randomseed=true",
                                                          "maximumiterations=5000","replicates=1",
                                                          "replicatetype=subsample", "-q", "betamultiplier=1"))
  List_Model[[i]] = MyModel
  TestAUC = c(TestAUC, MyModel@results["Test.AUC",1])
  MTSS_MyModel = c(MTSS_MyModel, MyModel@results["Maximum.test.sensitivity.plus.specificity.Cloglog.threshold",])
  MyPredict_Now = predict(MyModel, CurrentPredictors)
  values(MyPredict_Now)[values(MyPredict_Now) < MTSS_MyModel] = 0
  values(MyPredict_Now)[values(MyPredict_Now) > MTSS_MyModel] = 1
  List_Predict_Now[[i]] = MyPredict_Now
  MyPredict_Past = predict(MyModel, PastPredictors)
  values(MyPredict_Past)[values(MyPredict_Past) < MTSS_MyModel] = 0
  values(MyPredict_Past)[values(MyPredict_Past) > MTSS_MyModel] = 1
  List_Predict_Past[[i]] = MyPredict_Past
  
  if (as.numeric(MyModel@results["Test.AUC",1]) >= 0.8) {
    
    List_Predict_Now_AUC = append(List_Predict_Now_AUC, MyPredict_Now)
    List_Predict_Past_AUC = append(List_Predict_Past_AUC, MyPredict_Past)
    
  }
  
}

write.table(as.numeric(TestAUC), "1.1.1_AUC100Model.txt", quote = F, col.names = F, row.names = F)

png("1.1.2_BoxplotAUC100Model.png", height = 2100, width = 2100, res = 300)
# pdf("1.1.2_BoxplotAUC100Model.pdf", height = 7, width = 7)
boxplot(TestAUC, ylim = c(0,1), main = "AUC distribution", pch = 20, cex = 0.6)
dev.off()

Stack_Predict_Now = raster::stack(List_Predict_Now)
Stack_Predict_Past = raster::stack(List_Predict_Past)

# Remove past suitable pixels based on MESS and binary transformation based on the 90% threshold
Mean_Predict_Past = raster::stackApply(Stack_Predict_Past, indices = 100, fun = mean)
Mean_Predict_Past[Stack_MESS_Past_min$index_7 == 0] = 0
Mean_Predict_Past[Mean_Predict_Past < 0.9] = 0
Mean_Predict_Past[Mean_Predict_Past >= 0.9] = 1

# Same for current suitability but without removing MESS-based pixels
Mean_Predict_Now = raster::stackApply(Stack_Predict_Now, indices = 100, fun = mean)
Mean_Predict_Now[Mean_Predict_Now < 0.9] = 0
Mean_Predict_Now[Mean_Predict_Now >= 0.9] = 1

writeRaster(x = Mean_Predict_Past, "0_Rasters/Mean_Predict_Past_MESS")
writeRaster(x = Mean_Predict_Now, "0_Rasters/Mean_Predict_Now_MESS")

######################

# Only 5 models are retained because AUC equal or above 0.8

Stack_Predict_Past_AUC = raster::stack(List_Predict_Past_AUC)

Mean_Predict_Past_AUC = raster::stackApply(Stack_Predict_Past_AUC, indices = length(List_Predict_Past_AUC), fun = mean)
Mean_Predict_Past_AUC[Stack_MESS_Past_min$index_7 == 0] = 0
Mean_Predict_Past_AUC[Mean_Predict_Past_AUC < 0.9] = 0
Mean_Predict_Past_AUC[Mean_Predict_Past_AUC >= 0.9] = 1

Stack_Predict_Now_AUC = raster::stack(List_Predict_Now_AUC)

Mean_Predict_Now_AUC = raster::stackApply(Stack_Predict_Now_AUC, indices = length(List_Predict_Now_AUC), fun = mean)
Mean_Predict_Now_AUC[Mean_Predict_Now_AUC < 0.9] = 0
Mean_Predict_Now_AUC[Mean_Predict_Now_AUC >= 0.9] = 1

writeRaster(x = Mean_Predict_Past_AUC, "0_Rasters/Mean_Predict_Past_MESS_AUC")
writeRaster(x = Mean_Predict_Now_AUC, "0_Rasters/Mean_Predict_Now_MESS_AUC")
