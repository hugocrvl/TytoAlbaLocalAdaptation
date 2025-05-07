
###############################

########### LOCALLY ########### 

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('../0_Functions.R', chdir = T)

###############################

# I need to specify which window each SNP belongs to

###############################

RDA = fread("", header = T) # RDAAdapt output from ./0_RDA/2_RDA.R
MAF = fread("", header = F)$V1 # MAF value for every SNP, computed with ./1_MAF.sh
RDA$MAF = MAF
coordinates.gs = read.table("", header = T)[, c(1:6)] # FST zscores to have the window coordinates 

temp = data.frame()
for (s in unique(coordinates.gs$scaffold)) {
  
  subset.coordinates.gs = coordinates.gs[coordinates.gs$scaffold == s, ]
  subset.coordinates.gs$modulo = seq(1, nrow(subset.coordinates.gs)) %% 5

  subset.coordinates.gs$ID.wza = paste0(seq(1, nrow(subset.coordinates.gs)), "_",
                                            gsub("Super-Scaffold_", "", subset.coordinates.gs$scaffold))
  
  temp = rbind(temp, subset.coordinates.gs)
   
}

coordinates.gs = temp
remove(temp)
remove(subset.coordinates.gs)

get.modulo.subset = function(dataframe.RDA, dataframe.modulo, modulo) {
  
  print(paste("modulo :", modulo))
  subset.modulo = dataframe.modulo[dataframe.modulo$modulo == modulo, ]
  output = data.frame()
  
  for (i in 1:nrow(subset.modulo)) {
    
    if (i %in% c(1, seq(5e2,1e4,5e2))) {print(paste0('[', i, '/', nrow(subset.modulo), ']'))}
    
    row = subset.modulo[i, ]
    s = row$scaffold
    subset.RDA = dataframe.RDA[dataframe.RDA$scaffold == s, ]
    index.snps = which(subset.RDA$position >= row$start & subset.RDA$position <= row$end)
    
    if (length(index.snps) != 0) {
      
      subset.RDA.window = subset.RDA[index.snps]
      temp = data.frame(pvalue = subset.RDA.window$p.values, ID = row$ID.wza, 
                        MAF = subset.RDA.window$MAF)
      output = rbind(output, temp)
      
    }
  }
  
  return(output)
  
}

WZA.input.1 = get.modulo.subset(RDA, coordinates.gs, modulo = 1)
fwrite(WZA.input.1, "", quote = F)

WZA.input.2 = get.modulo.subset(RDA, coordinates.gs, modulo = 2)
fwrite(WZA.input.2, "", quote = F)

WZA.input.3 = get.modulo.subset(RDA, coordinates.gs, modulo = 3)
fwrite(WZA.input.3, "", quote = F)

WZA.input.4 = get.modulo.subset(RDA, coordinates.gs, modulo = 4)
fwrite(WZA.input.4, "", quote = F)

WZA.input.0 = get.modulo.subset(RDA, coordinates.gs, modulo = 0)
fwrite(WZA.input.0, "", quote = F)
