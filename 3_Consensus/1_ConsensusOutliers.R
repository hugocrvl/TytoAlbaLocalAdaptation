
###############################

########### LOCALLY ########### 

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("P:/Biologie/PhD/0_Functions.R", chdir = T)

###############################

tmp = read.table("", header = T) # FST outliers from ../1_DosageScans/4_GetOutliers.R
coordinates = tmp[, c(1:6)]
FST.outliers = tmp[, -c(1:6)]
FST.outliers.common = rowSums(FST.outliers)
FST.outliers.common[FST.outliers.common >= 1] = 1

WZA = read.table("", header = T) # From ../2_WZA/4_MergingWindows.R
WZA.outliers = read.table("", header = F)$V1 # From ../2_WZA/4_MergingWindows.R

consensus.outliers = FST.outliers.common + WZA.outliers
table(consensus.outliers)

write.table(consensus.outliers, "", quote = F, row.names = F, col.names = F)

###############################

consensus.outliers[consensus.outliers != 2] = 0
FST.outliers[-which(consensus.outliers == 2), ] = 0
FST.outliers.consensus = FST.outliers

write.table(FST.outliers.consensus, "", quote = F,
            col.names = T, row.names = F)

###############################

# Converting window coordinates into bed format to detect consecutive outlier regions

consensus.outliers = read.table("", header = F) # Same data.fram as in line 22
coordinates = read.table("", header = T)[, c(1:6)] # FST outliers from ../1_DosageScans/4_GetOutliers.R

coordinates.consensus.outliers = coordinates[consensus.outliers == 2, 1:3]

coordinates.consensus.outliers$scaffold = gsub("Super-Scaffold_", "chr", coordinates.consensus.outliers$scaffold)
coordinates.consensus.outliers$start = as.integer(format(coordinates.consensus.outliers$start, scientific = FALSE))
coordinates.consensus.outliers$end = as.integer(format(coordinates.consensus.outliers$end, scientific = FALSE))

write.table(coordinates.consensus.outliers, "", quote = F, col.names = F, row.names = F, sep = "\t")

