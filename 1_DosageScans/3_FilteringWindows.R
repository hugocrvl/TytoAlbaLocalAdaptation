
###############################

########### LOCALLY ###########

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../0_Functions.R", chdir = T)

###############################

# Comparison between MAF and without MAF statistics
# Need to remove extra windows from the without MAF statistics
# in order to compare the two statistics (with and without MAF)

# Again, Tajima's D wasn't used in the manuscript.

###############################

FST.MAF = read.table("", header = T)
FST.woMAF = read.table("", header = T)

tD.MAF = read.table("", header = T)
tD.woMAF = read.table("", header = T)

# Get in which scaffold there are differences
table.scaffold = table(tD.MAF$scaffold) == table(tD.woMAF$scaffold)
unequal.nw = names(which(table.scaffold == FALSE))

# Remove extra lines to be able to compare the two
output.FST = data.frame()
output.tD = data.frame()
for (s in unique(tD.woMAF$scaffold)) {
  
  print(s)
  
  if (s %in% unequal.nw) {
    
    n.diff = nrow(tD.woMAF[tD.woMAF$scaffold == s,]) - nrow(tD.MAF[tD.MAF$scaffold == s,])
    upper.limit = nrow(tD.woMAF[tD.woMAF$scaffold == s, ]) - n.diff
    
    temp.tD = tD.woMAF[tD.woMAF$scaffold == s, ]
    temp.tD = temp.tD[c(1:upper.limit), ]
    output.tD = rbind(output.tD, temp.tD)
    
    temp.FST = FST.woMAF[FST.woMAF$scaffold == s, ]
    temp.FST = temp.FST[c(1:upper.limit), ]
    output.FST = rbind(output.FST, temp.FST)
    
  } else {
    
    temp.tD = tD.woMAF[tD.woMAF$scaffold == s, ]
    output.tD = rbind(output.tD, temp.tD)
    
    temp.FST = FST.woMAF[FST.woMAF$scaffold == s, ]
    output.FST = rbind(output.FST, temp.FST)
    
  }
  
  remove(temp.tD)
  remove(temp.FST)
  
}

FST.woMAF = output.FST
remove(output.FST)
tD.woMAF = output.tD
remove(output.tD)

###############################

# I can now compare the with and without MAF statistics to know which one I use

###############################

FST.corr.matrix = cor(as.matrix(FST.MAF[, -c(1:6)]), as.matrix(FST.woMAF[, -c(1:6)]), 
                      use = "complete.obs")
mean(diag(FST.corr.matrix))

###############################

# I can use the woMAF statistics as they are strongly correlated with the MAF ones
# Now I need to remove the windows with too few SNPs to have a certain level of info

###############################

FST.woMAF = read.table("", header = T)
tD.woMAF = read.table("", header = T)

m.n.snps = mean(FST.woMAF$n.snps, na.rm = TRUE)
sd.n.snps = sd(FST.woMAF$n.snps, na.rm = TRUE)
threshold = ceiling(m.n.snps - (2*sd.n.snps))
index.window.underT = which(FST.woMAF$n.snps <= threshold)
# length(index.window.underT) = 2524

# I remove windows with less than 250 SNPs (being 2SD below the mean) per window
FST.woMAFremovedNA = FST.woMAF[-index.window.underT, ]
tD.woMAFremovedNA = tD.woMAF[-index.window.underT, ]

# At the end, I use a data set of 52429 windows
write.table(FST.woMAFremovedNA, "", quote = F, col.names = T, row.names = F)
write.table(tD.woMAFremovedNA, "", quote = F, col.names = T, row.names = F)
