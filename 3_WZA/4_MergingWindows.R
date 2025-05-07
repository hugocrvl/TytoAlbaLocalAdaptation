
###############################

########### LOCALLY ########### 

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('../0_Functions.R', chdir = T)

###############################

# I'll merge the output from the 5 different WZA and sort them

###############################

WZA.output.merged = read.csv("", header = T) # Merged WZA output from ./3_WZA.sh
coordinates.gs = read.table("", header = T)[, c(1:6)] # FST zscores to have the window coordinates 

list = strsplit(WZA.output.merged$gene, '_')
scaffold.vector = c()
window.vector = c()

for (i in 1:length(list)) {
  
  if (i %in% c(seq(0, 5e4, 1e4))) {print(i)}
  scaffold.vector = c(scaffold.vector, list[i][[1]][2])
  window.vector = c(window.vector, list[i][[1]][1])
  
}

WZA.output.merged$scaffold = as.numeric(scaffold.vector)
WZA.output.merged$window = as.numeric(window.vector)
WZA.output.merged = WZA.output.merged[order(WZA.output.merged$scaffold, decreasing = F), ]

WZA = data.frame()
for (s in unique(WZA.output.merged$scaffold)) {
  
  subset.WZA.output.merged = WZA.output.merged[WZA.output.merged$scaffold == s, ]
  WZA = rbind(WZA, subset.WZA.output.merged[order(subset.WZA.output.merged$window, decreasing = F), ])
  
}

write.table(WZA, "", quote = F, col.names = T, row.names = F)

###############################

# To get WZA outliers

WZA = read.table("", header = T) # Output from line 44
WZA.outliers = rep(0, nrow(WZA))
threshold = mean(-log10(WZA$Z_pVal)) + (2*sd(-log10(WZA$Z_pVal)))
WZA.outliers[which(-log10(WZA$Z_pVal) >= threshold)] = 1

write.table(WZA.outliers, "", quote = F, col.names = F, row.names = F)
