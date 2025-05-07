
############################################

########### NOT TO RUN - CLUSTER ###########

############################################

setwd("")
source("../../0_Functions.R", chdir = T)

library("vegan")
library("robust")
library("qvalue")

###############################

print("Opening of the VCF")
VCF = read.vcf("", convert.chr = F)
dosage = as.matrix(VCF)
print("Done")
explanatory = read.table("", header = T) # Explanatory variables from ./1_Explanatory.R

print("Starting RDA")
rda = rda(dosage ~ bio2 + bio6 + bio7 + bio8 + bio15 + bio17 + bio19, data = explanatory)
RsquareAdj(rda)$adj.r.squared
print("Done")

rdadapt = function(rda, K) {
  
  zscores =rda$CCA$v[,1:as.numeric(K)]
  resscale = apply(zscores, 2, scale)
  resmaha = covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda = median(resmaha)/qchisq(0.5,df=K)
  reschi2test = pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval = qvalue(reschi2test)
  q.values_rdadapt = qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))

}

result.rdadapt = rdadapt(rda, K = 5)
result.rdadapt = data.frame(scaffold = VCF@snps$chr, position = VCF@snps$pos, result.rdadapt)
write.table(result.rdadapt, "", quote = F, col.names = T, row.names = F)

###############################

# Permutation test for model significance
anova(rda, permutations = how(nperm = 999))

# Permutation test for axis significance
anova(rda, by = "axis", permutations = how(nperm = 100))
