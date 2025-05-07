
###############################

########### LOCALLY ########### 

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("P:/Biologie/PhD/0_Functions.R", chdir = T)

###############################

library("SNPRelate")

sample.info = read.table("", header = T) # Individual and population IDs
GDS = snpgdsOpen("") # Same dataset as with FST scans
ID.SNP = snpgdsSNPList(GDS, sample.id = NULL)
populations = substr(sample.info$ID, 1, 2)

list.colour.bg = colour.population[as.factor(populations)]
list.pch = c(21:25, 21:24)[as.factor(populations)]

###############################

# PCA on the inversion

###############################

ID.SNP.sc22 = ID.SNP[ID.SNP$chromosome == "Super-Scaffold_22", ]
threshold = 14e6

PCA.inversion = snpgdsPCA(GDS, autosome.only = FALSE, verbose = TRUE, eigen.cnt = 0, remove.monosnp = FALSE,
                snp.id = ID.SNP.sc22$snp.id[ID.SNP.sc22$position <= threshold])

###############################

# PCA on the whole genome

###############################

PCA.wholeG = snpgdsPCA(GDS, autosome.only = FALSE, verbose = TRUE, eigen.cnt = 0, remove.monosnp = FALSE)
