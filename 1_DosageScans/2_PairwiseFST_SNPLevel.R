
############################################

setwd("")
library(hierfstat)
library(gaston)

pathVCF = ""

print("Opening VCF")
VCF = read.vcf(pathVCF, convert.chr = F)
dosage = as.matrix(VCF)
print("Done")

populations = substr(VCF@ped$id, 1, 2)
# All possible population pairs of your VCF, our case:
column.names = c("AE.CH", "AE.DK", "AE.FR", "AE.GR", "AE.IS",
                 "AE.IT", "AE.PT", "AE.SB", "CH.DK", "CH.FR",
                 "CH.GR", "CH.IS", "CH.IT", "CH.PT", "CH.SB",
                 "DK.FR", "DK.GR", "DK.IS", "DK.IT", "DK.PT",
                 "DK.SB", "FR.GR", "FR.IS", "FR.IT", "FR.PT",
                 "FR.SB", "GR.IS", "GR.IT", "GR.PT", "GR.SB",
                 "IS.IT", "IS.PT", "IS.SB", "IT.PT", "IT.SB",
                 "PT.SB")
pairwise.SNP = matrix(-9, nrow = ncol(dosage), ncol = length(column.names))
colnames(pairwise.SNP) = column.names

for (i in 1:ncol(dosage)) {
  
  if (i %in% c(1, seq(5e5, 12e6, 5e5))) {print(i)}
  
  column = dosage[, i]
  fs = fs.dosage(column, populations)
  pairwise = fs$Fst2x2
  
  for (j in 1:length(column.names)) {
    
    comparison = column.names[j]
    pop1 = substr(comparison, 1, 2)
    pop2 = substr(comparison, 4, 5)
    
    pairwise.SNP[i, comparison] = pairwise[pop1, pop2]
    
  }
  
}

write.table(pairwise.SNP, "", quote = F, col.names = T, row.names = F)