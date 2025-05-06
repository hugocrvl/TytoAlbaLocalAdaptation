
###############################

########### LOCALLY ########### 

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library("hierfstat")
library("gaston")

# Function to compute Tajima's D, Pi or population-specific FST along the genome
# It uses overlapping sliding windows of 100kb, with a step of 20kb,
# but these parameters can be changed

dosage.function = function(pathVCF, window.size, window.step, index = c("Pi", "FST", "tajimaD")) {
  
  if (index == "Pi") {print("Pi will be computed")} else if (index == "FST") 
    {print("Population-specific FST will be computed")} else {print("Tajima's D will be computed")}
  
  print("Opening VCF")
  VCF = read.vcf(pathVCF, convert.chr = F)
  dosage = as.matrix(VCF)
  print("Done")
  
  populations = substr(VCF@ped$id, 1, 2)
  positions = VCF@snps$pos
  scaffolds = VCF@snps$chr
  output = data.frame()
  
  for (scaffold in unique(scaffolds)) {
    
    print(paste0(scaffold, " [", which(unique(scaffolds) == scaffold), 
                 "/", length(unique(scaffolds)), "]"))
    subset.dosage = dosage[, VCF@snps$chr == scaffold]
    subset.positions = positions[VCF@snps$chr == scaffold]
    highest.position = tail(subset.positions, n=1)
    first = TRUE
    
    if (index == "FST") {
      
      start = 0
      end = start + window.size
      start.vector = c()
      end.vector = c()
      n.snps.vector = c()
      first.SNP.vector = c()
      last.SNP.vector = c()
      
      stats.matrix = matrix(-999, ncol = sum(lower.tri(matrix(0, ncol = length(unique(populations)),
                                                        nrow = length(unique(populations))), diag = T)))
      colnames(stats.matrix) = colnames = c(unique(populations),
                                            "AE.CH", "AE.DK", "AE.FR", "AE.GR", "AE.IS",
                                            "AE.IT", "AE.PT", "AE.SB", "CH.DK", "CH.FR",
                                            "CH.GR", "CH.IS", "CH.IT", "CH.PT", "CH.SB",
                                            "DK.FR", "DK.GR", "DK.IS", "DK.IT", "DK.PT",
                                            "DK.SB", "FR.GR", "FR.IS", "FR.IT", "FR.PT",
                                            "FR.SB", "GR.IS", "GR.IT", "GR.PT", "GR.SB",
                                            "IS.IT", "IS.PT", "IS.SB", "IT.PT", "IT.SB",
                                            "PT.SB")
      
      continue = TRUE
      
      while(continue) {
        
        index.SNP.window = which(subset.positions >= start & subset.positions <= end)
        n.snps.vector = c(n.snps.vector, length(index.SNP.window))
        
        if (length(index.SNP.window) == 0) {
          
          start.vector = c(start.vector, start)
          end.vector = c(end.vector, end)
          
          first.SNP.vector = c(first.SNP.vector, NA)
          last.SNP.vector = c(last.SNP.vector, NA)
          
          stats.matrix = rbind(stats.matrix, rep(NA, n = sum(lower.tri(matrix(0, ncol = length(unique(populations)),
                                                                              nrow = length(unique(populations))), diag = T))))
          
        } else {
        
          start.vector = c(start.vector, start)
          end.vector = c(end.vector, end)
          
          window = subset.dosage[, index.SNP.window]
          subset.positions.window = subset.positions[index.SNP.window]
          
          first.SNP = as.numeric(subset.positions.window[1])
          first.SNP.vector = c(first.SNP.vector, first.SNP)
          
          last.SNP = as.numeric(tail(subset.positions.window, n=1))
          last.SNP.vector = c(last.SNP.vector, last.SNP)
          
          fs = fs.dosage(window, pop = populations)
          popspe.vector = diag(fs$FsM)
          pairwise.vector = fs$FsM[lower.tri(fs$FsM, diag = F)]

          stats.matrix = rbind(stats.matrix, c(popspe.vector, pairwise.vector))

        }
        
        start = start + window.step
        end = start + window.size
        if (end >= highest.position) {break}
        
      }
      
      temp = data.frame(scaffold = scaffold, start = start.vector, end = end.vector,
                        n.snps = n.snps.vector, first.SNP = first.SNP.vector,
                        last.SNP = last.SNP.vector, stats.matrix[-1, ])
      
    }
    
    if (index != "FST") {
      
      for (population in unique(populations)) {
        
        print(population)
        subset.dosage.ID = subset.dosage[grepl(population, rownames(subset.dosage)), ]
        
        start = 0
        end = start + window.size
        start.vector = c()
        end.vector = c()
        n.snps.vector = c()
        first.SNP.vector = c()
        last.SNP.vector = c()
        stats.vector = c()
        
        continue = TRUE
        
        while(continue) {
          
          index.SNP.window = which(subset.positions >= start & subset.positions <= end)
          n.snps.vector = c(n.snps.vector, length(index.SNP.window))
          
          if (length(index.SNP.window) == 0) {
            
            start.vector = c(start.vector, start)
            end.vector = c(end.vector, end)
            
            first.SNP.vector = c(first.SNP.vector, NA)
            last.SNP.vector = c(last.SNP.vector, NA)
            
            stats.vector = c(stats.vector, NA)
            
          } else {
            
            window = subset.dosage.ID[, index.SNP.window]
            subset.positions.window = subset.positions[index.SNP.window]
            
            start.vector = c(start.vector, start)
            end.vector = c(end.vector, end)
            
            first.SNP = as.numeric(subset.positions.window[1])
            first.SNP.vector = c(first.SNP.vector, first.SNP)
            
            last.SNP = as.numeric(tail(subset.positions.window, n=1))
            last.SNP.vector = c(last.SNP.vector, last.SNP)
            
            if (index == "Pi") {stats.vector = c(stats.vector, pi.dosage(window, L = last.SNP - first.SNP))} else {
              stats.vector = c(stats.vector, TajimaD.dosage(window))}
            
          }
          
          start = start + window.step
          end = start + window.size
          if (end >= highest.position) {break}
          
        }
        
        if (first) {
          temp = data.frame(scaffold = scaffold, start = start.vector, end = end.vector,
                            n.snps = n.snps.vector, first.SNP = first.SNP.vector,
                            last.SNP = last.SNP.vector, stats.vector = stats.vector)
          colnames(temp) = c("scaffold", "start", "end", "n.snps", "first.SNP", 
                             "last.SNP", paste0(population))
          first = FALSE
        } else {
          
          temp = cbind(temp, stats.vector)
          colnames(temp)[ncol(temp)] = paste0(population)
          
          }
        
        # temp = data.frame(scaffold = scaffold, start = start.vector, end = end.vector,
        #                   n.snps = n.snps.vector, first.SNP = first.SNP.vector,
        #                   last.SNP = last.SNP.vector)
        
        
      }
      
    }
    
    output = rbind(output, temp)
    
  }
  
  return(output)

}

# Now, let's use this function to compute population-specific FST, Tajima's D and Pi for VCF with or without MAF filtering

# Population-specific FST

FST.MAF = dosage.function(pathVCF = "C:/Users/Hugo Corval/Documents/Complete_EU_filters_masked.vcf.gz",
                      window.size = 100000, window.step = 20000, index = "FST")
write.table(FST.MAF, "1.1.1_FST.MAF.txt", quote = F, col.names = T, row.names = F)
FST.woMAF = dosage.function(pathVCF = "C:/Users/Hugo Corval/Documents/Complete_EU_filters_masked_woMAF.vcf.gz",
                            window.size = 100000, window.step = 20000, index = "FST")
write.table(FST.woMAF, "1.1.2_FST.woMAF.txt", quote = F, col.names = T, row.names = F)

# Tajima's D

tD.MAF = dosage.function(pathVCF = "C:/Users/hcorval/Documents/Complete_EU_filters_masked.vcf.gz",
                         window.size = 100000, window.step = 20000, index = "tajimaD")
write.table(tD.MAF, "1.2.1_tDMAF.txt", quote = F, col.names = T, row.names = F)
tD.woMAF = dosage.function(pathVCF = "C:/Users/hcorval/Documents/Complete_EU_filters_masked_woMAF.vcf.gz",
                           window.size = 100000, window.step = 20000, index = "tajimaD")
write.table(tD.woMAF, "1.2.2_tDwoMAF.txt", quote = F, col.names = T, row.names = F)

# Pi

# Pi.MAF = dosage.function(pathVCF = "/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/LocalAdaptation/Complete_EU_filters_masked.vcf",
#                          window.size = 100000, window.step = 20000, index = "Pi")
# write.table(Pi.MAF, "1.3.1_PiMAF.txt", quote = F, col.names = T, row.names = F)
# Pi.woMAF = dosage.function(pathVCF = "/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/LocalAdaptation/phasing/Complete_EU_filters_masked_woMAF.vcf",
#                            window.size = 100000, window.step = 20000, index = "Pi")
# write.table(Pi.woMAF, "1.3.2_PiwoMAF.txt", quote = F, col.names = T, row.names = F)

############################################

setwd("/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/LocalAdaptation/1_DosageScans/")
library(hierfstat)
library(gaston)

pathVCF = "/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/LocalAdaptation/phasing/Complete_EU_filters_masked_woMAF.vcf.gz"

print("Opening VCF")
VCF = read.vcf(pathVCF, convert.chr = F)
dosage = as.matrix(VCF)
print("Done")

populations = substr(VCF@ped$id, 1, 2)
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

write.table(pairwise.SNP, "1.3_PairwiseFST.SNP.txt", quote = F, col.names = T, row.names = F)
