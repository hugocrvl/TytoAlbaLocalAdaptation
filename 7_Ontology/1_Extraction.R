
###############################

########### LOCALLY ########### 

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../0_Functions.R", chdir = T)

###############################

# Extraction of the background list of genes for the Gene Enrichment Analysis 
# All the genes spanning FST/WZA windows

coordinates = read.table("", header = T)[, c(1:6)] # FST outliers from ../1_DosageScans/4_GetOutliers.R for windows coordinates
table.genes = read.table("", header = T) # Gene annotation from the annotated barn owl genome (GenBank Assembly Accession: GCA_018691265.1)

genes.background = c()

for (s in unique(coordinates$scaffold)) {
  
  print(s)
  table.genes.scaffold = table.genes[table.genes$scaffold == s, ]
  coordinates.scaffold = coordinates[coordinates$scaffold == s, ]
  
  for (i in 1:nrow(coordinates.scaffold)) {
    
    start = coordinates.scaffold$start[i]
    end = coordinates.scaffold$end[i]
    
    gene.start.within = which(table.genes.scaffold$start >= start &
                                table.genes.scaffold$start <= end)
    gene.end.within = which(table.genes.scaffold$end >= start &
                              table.genes.scaffold$end <= end)
    gene.bigger.window = which(table.genes.scaffold$start <= start &
                                 table.genes.scaffold$end >= end)
    all.index = unique(c(gene.start.within, gene.end.within, gene.bigger.window))
    
    if (length(all.index) != 0) {
      
      genes.background = c(genes.background, table.genes.scaffold$Gene[all.index])    
      
    }
    
  }
  
}

genes.background = unique(genes.background)

write.table(genes.background, "1.1_GenesBackground.txt", quote = F, col.names = F, row.names = F)

###############################

# Extraction of genes falling in outlier windows 
# First, per population
# Second, merged list of genes

coordinates = read.table("", header = T)[, c(1:6)] #  FST outliers from ../1_DosageScans/4_GetOutliers.R for windows coordinates
consensus.outliers = read.table("", header = F)$V1 # List of consensus outliers from ../1_Consensus/1_ConsensusOutliers.R
population.matrix.consensus.outliers = read.table("", header = T) # Matrixc of consensus outliers per population from ../1_Consensus/1_ConsensusOutliers.R
table.genes = read.table("", header = T) # Gene annotation from the annotated barn owl genome (GenBank Assembly Accession: GCA_018691265.1)

coordinates.outliers = coordinates[consensus.outliers == 2, ]

output = data.frame()

for (p in colnames(population.matrix.consensus.outliers)) {
  
  print(p)
  population.coordinates = coordinates[population.matrix.consensus.outliers[, p] == 1, ]
  
  population.genes = data.frame()
  
  for (i in 1:nrow(population.coordinates)) {
    
    row.coordinates = population.coordinates[i, ]
    genes.scaffold = table.genes[table.genes$scaffold == row.coordinates$scaffold, ]
    
    gene.start.within = which(genes.scaffold$start >= row.coordinates$start &
                                genes.scaffold$start <= row.coordinates$end)
    gene.end.within = which(genes.scaffold$end >= row.coordinates$start &
                              genes.scaffold$end <= row.coordinates$end)
    gene.bigger.window = which(genes.scaffold$start <= row.coordinates$start &
                            genes.scaffold$end >= row.coordinates$end)
    
    all.index = unique(c(gene.start.within, gene.end.within, gene.bigger.window))
    extracted.genes = genes.scaffold[all.index, ]
    
    if (length(all.index) != 0) {
      
      temp = data.frame(extracted.genes, population = p)
      population.genes = rbind(population.genes, temp)
    
    }
    
  }
  
  output = rbind(output, population.genes)
  
}

write.table(output, "", quote = F, row.names = F, col.names = T) # Table with all outlier genes per population
write.table(unique(output$Gene), "", quote = F, row.names = F, col.names = F) # Unique outlier genes across populations
