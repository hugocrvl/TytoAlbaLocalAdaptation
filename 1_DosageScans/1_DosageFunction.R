
###############################

########### LOCALLY ########### 

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library("hierfstat")
library("gaston")

# Function to compute Tajima's D, Pi or population-specific and population-pair FST along the genome
# It uses overlapping sliding windows of 100kb, with a step of 20kb,
# but these parameters can be changed.
# In the manuscript, Tajima's D and Pi weren't used.

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
        
      }
      
    }
    
    output = rbind(output, temp)
    
  }
  
  return(output)

}
