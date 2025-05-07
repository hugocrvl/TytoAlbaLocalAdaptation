
###############################

########### LOCALLY ########### 

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../0_Functions.R", chdir = T)

###############################

set.seed(1)

######### Load Required Libraries #########

library(hierfstat)   # Genetic diversity and population structure analyses
library(SNPRelate)   # Functions for handling SNP data
library(ade4)        # Multivariate analysis (e.g., PCA)
library(pegas)       # Analysis of population genetics
library(gplots)      # Enhanced plotting functions
library(randomcoloR) # Generate random colors
library(viridis)   # Color palettes for visualization

COLOR = randomColor(13)

##### Load  Data #####

DATA = snpgdsOpen("")  # Open SNP dataset, same as for ../5_PCA/1_MakingPCA.R
gen_mat = snpgdsGetGeno("") # Extract genotype matrix, same gds as previous line
dim(gen_mat) # Check dimensions

##### Genetic Diversity Analysis #####

# Inbreeding Coefficient Calculation
ImbCoef = snpgdsIndInb(DATA, autosome.only = FALSE)
IndId = ImbCoef$sample.id  # Individual IDs
Pop = substr(ImbCoef$sample.id, 1, 2)  # Extract population codes

# Heterozygosity Calculation
IndHet = rowSums(gen_mat == 1, na.rm = TRUE) / rowSums(!is.na(gen_mat))
boxplot(IndHet ~ Pop)  # Boxplot of heterozygosity per population
plot(IndHet[order(Pop, IndHet)], pch = 20, col = COLOR[as.factor(Pop)])

## Population Scale Analysis

# We need to downsample populations to have an homogeneous sampling of alleles
# We did not downsampled France : only 4 samples, too few individuals, diversity similar to CH or DK
# We have to reduce to 5 to match SB sampling size

SampleInd=function(PopVector, PopId, NombreInd){
  sample(c(1:length(PopVector))[as.factor(PopVector)==PopId], size = NombreInd)
}

# Define populations of interest
Pop2Keep = c("AE","CH","DK","GR","IS","IT","PT","SB")
Polymorphic = NULL  # To store counts of polymorphic sites
Private = NULL  # To store counts of private alleles
Rare = NULL  # To store counts of rare alleles

gc() # Free memory

for (Rep in 1:10){  # Repeat sampling process 10 times to ensure robustness
  cat("\n")
  cat(paste0("Bootstrap ", Rep, "\n"))
  Ind2Keep = c()
  PopMatrix = list()
  
  # Downsample populations to ensure equal representation
  for (p in unique(Pop2Keep)){
    SelectedInd = SampleInd(Pop, p, 5)  # Sample 5 individuals per population
    Ind2Keep = c(Ind2Keep, SelectedInd)
    PopMatrix[[p]] = gen_mat[SelectedInd,]
  }
  
  # Include France (FR) even with fewer individuals
  p = "FR"
  SelectedInd = which(Pop == p)
  Ind2Keep = c(Ind2Keep, SelectedInd)
  PopMatrix[[p]] = gen_mat[SelectedInd,]
  
  # Create a full genotype matrix
  FullMat = do.call(rbind, PopMatrix)

  # Identify polymorphic sites per population
  IsPolymorph = NULL # Each column will correspond to a population, each row to a SNP, filled with logical value (TRUE for polymorphic, FALSE for monomorphic)
  cat("\n")
  cat("Find polymorphic sites per population\n")
  for (p in c(Pop2Keep, "FR")){
    print(p)
    IsPolymorph = cbind(IsPolymorph, apply(PopMatrix[[p]], 2, function(x){length(unique(x[!is.na(x)])) > 1}))
  }
  cat("\n")
  
  # Calculate allele frequencies
  cat("Computing allelic frequencies\n")
  Freq = c()  # Initialize an empty vector to store frequencies
  chunk_size = 1e6  # Process 1 million SNPs at a time (adjust as needed)
  num_chunks = ceiling(ncol(FullMat) / chunk_size)  # Calculate how many chunks are needed
  
  for (i in 1:num_chunks) {
    start_idx = (i - 1) * chunk_size + 1
    end_idx = min(i * chunk_size, ncol(FullMat))  # Avoid exceeding matrix dimensions
    
    cat("Processing SNPs from", start_idx, "to", end_idx, "\n")  # Track progress
    
    FullMat_chunk = FullMat[, start_idx:end_idx]  # Extract the chunk
    Freq_chunk = colSums(FullMat_chunk, na.rm = TRUE) / (colSums(!is.na(FullMat_chunk)) * 2)  # Compute allele frequency
    
    Freq = c(Freq, Freq_chunk)  # Store results
    
    rm(FullMat_chunk)  # Free memory
    gc()  # Trigger garbage collection
  }
  
  IsRare = (Freq < 0.05 | Freq > 0.95)  # Define rare alleles as minor allele frequency < 5%
  
  # Compute diversity statistics for each population
  bootstrap.polymorphic = c()
  bootstrap.private = c()
  bootstrap.rare = c()
  for (i in 1:9){
    bootstrap.polymorphic = c(bootstrap.polymorphic, sum(IsPolymorph[, i] == 1))  # Polymorphic sites
    bootstrap.private = c(bootstrap.private, sum(IsPolymorph[, i] == 1 & rowSums(IsPolymorph[, -i]) == 0))  # Private alleles
    bootstrap.rare = c(bootstrap.rare, sum(IsPolymorph[, i] == 1 & IsRare == TRUE))  # Rare alleles
  }
  
  Polymorphic = rbind(Polymorphic, bootstrap.polymorphic)
  Private = rbind(Private, bootstrap.private)
  Rare = rbind(Rare, bootstrap.rare)
  
  gc()  # Trigger garbage collection
  
}

# Assign column names to the diversity matrices, ensuring consistency across datasets
colnames(Polymorphic) = colnames(Private) = colnames(Rare) = c(Pop2Keep, "FR")

# Save the diversity statistics to text files for further analysis
write.table(Polymorphic, "1.1_Polymorphic.txt", quote = F, row.names = F, col.names = T)
write.table(Private, "1.2_Private.txt", quote = F, row.names = F, col.names = T)
write.table(Rare, "1.3_Rare.txt", quote = F, row.names = F, col.names = T)

###############################

Polymorphic = read.table("1.1_Polymorphic.txt", header = T)
Private = read.table("1.2_Private.txt", header = T)
Rare = read.table("1.3_Rare.txt", header = T)
Pop2Keep = c("AE","CH","DK","GR","IS","IT","PT","SB")

# Compute mean and standard deviation for each diversity metric per population
mean.polymorphic = colMeans(Polymorphic)
sd.polymorphic = round(apply(Polymorphic, 2, sd, na.rm = TRUE))

mean.private = colMeans(Private)
sd.private = round(apply(Private, 2, sd, na.rm = TRUE))

mean.rare = colMeans(Rare)
sd.rare = round(apply(Rare, 2, sd, na.rm = TRUE))

DiversityTable = data.frame(
  Population = c(Pop2Keep, "FR"),
  SampleSize = as.numeric(table(Pop)[c(Pop2Keep, "FR")]),
  Mean_Polymorphic = round(colMeans(Polymorphic)),
  SD_Polymorphic = round(apply(Polymorphic, 2, sd, na.rm = TRUE)),
  Mean_Private = round(colMeans(Private)),
  SD_Private = round(apply(Private, 2, sd, na.rm = TRUE)),
  Mean_Rare = round(colMeans(Rare)),
  SD_Rare = round(apply(Rare, 2, sd, na.rm = TRUE))
)

write.table(DiversityTable, "1.5_DiversityTable.txt", quote = F, row.names = T, col.names = T)

###############################

#### TO RUN ON THE CLUSTER ####

###############################

#### Population-specific FST #####

setwd("")

library(hierfstat)  # Used for computing FST and other genetic diversity metrics
library(gaston)     # Used for handling VCF files and genotype dosage calculations

# Read the VCF file containing filtered genetic data, same as for the diversity table
VCF = read.VCF("", convert.chr = F)
dosage = as.matrix(VCF)

# Extract population identifiers from individual IDs in the VCF file
populations = substr(VCF@ped$id, 1, 2)

print("Computing whole genome population-specific FST")
# Compute population-specific FST values using genotype dosage data
# fs.dosage calculates FST for each population based on the genetic variance
FST = fs.dosage(dosage, populations)
FST.matrix = FST$FsM

# overall FST can be retrieved by FST$Fs["Fst", "All"]

# Save the computed FST matrix to a text file for further analysis
write.table(FST.matrix, "1.6_FST.txt", quote = F, row.names = F, col.names = T)

###############################

########### LOCALLY ########### 

###############################

DiversityTable = read.table("1.5_DiversityTable.txt", header = T)
FST = as.matrix(read.table("1.6_FST.txt", header = T))
rownames(FST) = colnames(FST) 

populations.in.order = c("AE","CH","DK","GR","IS","IT","PT","SB", "FR")
FST = FST[populations.in.order, populations.in.order]
population.specific.FST = diag(FST)

TableManuscript = cbind(DiversityTable, WholeG_PopSpecificFST = round(population.specific.FST, 3))
TableManuscript.PopulationOrder = TableManuscript[order(TableManuscript$Population), ]

write.table(TableManuscript.PopulationOrder, "1.7_Diversity_TableManuscript.txt", quote = F, col.names = T)
