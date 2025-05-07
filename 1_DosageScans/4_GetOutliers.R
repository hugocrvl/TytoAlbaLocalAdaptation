
###############################

########### LOCALLY ###########

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../0_Functions.R", chdir = T)

###############################

# Script to get the outlier windows for FST.
# Again, Tajima's wasn't used in the manuscript.

# I compute the z-scores for each population and save into files

###############################

FST = read.table("", header = T) # output from ./3_FilteringWindows.R
tD = read.table("", header = T) 

get.zscores = function(dataframe) {
  
  temp = dataframe[, c(1:6)]
  columns = colnames(dataframe[, -c(1:6)])
  
  for (c in columns) {
    stats = dataframe[, c]
    m = mean(stats)
    s = sd(stats)
    z = (stats-m) / s
    temp = cbind(temp, z)
  }
  
  colnames(temp)[-c(1:6)] = columns
  return(temp)
  
}

FST.zscores = get.zscores(FST)
tD.zscores = get.zscores(tD)
write.table(FST.zscores, "", quote = F, col.names = T, row.names = F)
write.table(tD.zscores, "", quote = F, col.names = T, row.names = F)

###############################

# For every population, I detect outliers 

###############################

FST.zscores = read.table("", header = T) # z-scores FST computed above
tD.zscores = read.table("", header = T)

get.outliers.index = function(dataframe, side = c("right", "left"), figure = c(TRUE, FALSE)) {
  
  zscores.vector = c()
  
  for (p in colnames(dataframe[, c(7:15)])) {zscores.vector = c(zscores.vector, dataframe[, p])}
  
  m = mean(zscores.vector)
  s = sd(zscores.vector)
  print(m);print(s)
  if (side == "right") {threshold = m + (2*s)} else {threshold = m - (2*s)}
  print(threshold)
  if (figure) {
    
    hist(zscores.vector, breaks = 100)
    abline(v = threshold, col = "red", lty = 2, lwd = 2)
    
  }
  
  temp = dataframe[, c(1:6)]
  for (p in colnames(dataframe[, c(7:15)])) {
    
    binary = rep(0, nrow(dataframe))
    if (side == "right") {binary[dataframe[, p] >= threshold] = 1} else {
      binary[dataframe[, p] <= threshold] = 1}
    
    temp = cbind(temp, binary)
    colnames(temp)[ncol(temp)] = paste(p)
    
  }
  
  return(temp)
  
}

FST.outliers = get.outliers.index(FST.zscores, side = "right", figure = FALSE)
tD.outliers = get.outliers.index(tD.zscores, side = "left", figure = FALSE)
write.table(FST.outliers, "", quote = F, col.names = T, row.names = F)
write.table(tD.outliers, "", quote = F, col.names = T, row.names = F)
