
###############################

########### LOCALLY ###########

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../0_Functions.R", chdir = T)

###############################

# I need to remove extra windows from the without MAF statistics
# in order to compare the two statistics (with and without MAF)

###############################

FST.MAF = read.table("1.1.1_FST.MAF.txt", header = T)
FST.woMAF = read.table("1.1.2_FST.woMAF.txt", header = T)

tD.MAF = read.table("1.2.1_tDMAF.txt", header = T)
tD.woMAF = read.table("1.2.2_tDwoMAF.txt", header = T)

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
# 0.8797208

png("../SupplementaryMaterial/FSTRelationFSTwoMAF.png", height = 3000, width = 3000, res = 300)

plot(as.matrix(FST.MAF[, -c(1:6)]), as.matrix(FST.woMAF[, -c(1:6)]), xlab = "Population-specific FST with MAF",
     ylab = "Population-specific FST without MAF", col = add.alpha("black", alpha = 0.80), pch = 20)

dev.off()

ctD.corr.matrix = cor(as.matrix(tD.MAF[, -c(1:6)]), as.matrix(tD.woMAF[, -c(1:6)]), 
                     use = "complete.obs")
mean(diag(tD.corr.matrix))
# 0.8707618
plot(as.matrix(tD.MAF[, -c(1:6)]), as.matrix(tD.woMAF[, -c(1:6)]), xlab = "Tajima's D with MAF",
     ylab = "Tajima's D without MAF", col = add.alpha("black", alpha = 0.80), pch = 20)

###############################

# I can use the woMAF statistics as they are strongly correlated with the MAF ones
# Now I need to remove the windows with too few SNPs to have a certain level of info

###############################

FST.woMAF = read.table("1.1.2_FstwoMAF.txt", header = T)
tD.woMAF = read.table("1.2.2_tDwoMAF.txt", header = T)

# cor(tD.woMAF$n.snps, FST.woMAF$n.snps, use = "complete.obs")
# cor = 1, does not matter from which data set I extract the number of SNP
m.n.snps = mean(FST.woMAF$n.snps, na.rm = TRUE)
sd.n.snps = sd(FST.woMAF$n.snps, na.rm = TRUE)
hist(FST.woMAF$n.snps, breaks = 100)
abline(v = m.n.snps, col = "red", lty = 2, lwd = 2)
abline(v = c(m.n.snps - (2*sd.n.snps)), col = "black", lty = 2, lwd = 2)

threshold = ceiling(m.n.snps - (2*sd.n.snps))
# 250
index.window.underT = which(FST.woMAF$n.snps <= threshold)
# length(index.window.underT) = 2524

# I remove windows with less than 250 SNPs (being 2SD below the mean) per window
FST.woMAFremovedNA = FST.woMAF[-index.window.underT, ]
tD.woMAFremovedNA = tD.woMAF[-index.window.underT, ]

# At the end, I use a data set of 52429 windows
write.table(FST.woMAFremovedNA, "2.1_FSTwoMAFremovedNA.txt", quote = F, col.names = T, row.names = F)
write.table(tD.woMAFremovedNA, "2.2_tDwoMAFremovedNA.txt", quote = F, col.names = T, row.names = F)
