
###############################

########### LOCALLY ###########

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../0_Functions.R", chdir = T)

###############################

# Figures of population-specific FST for every population

###############################

FST.zscores = read.table("3.1_FST.zscores.txt", header = T)
FST.outliers = read.table("3.3.1_FST.outliers.txt", header = T)

middle = get.middle.position(FST.zscores$start, FST.zscores$end)
absolute.position = get.absolute.position(middle, FST.zscores$scaffold)
col.scaffold = get.scaffold.colour(FST.zscores$scaffold)
thick.mark = get.thick.mark(middle, FST.zscores$scaffold)

for (p in colnames(FST.outliers)[-c(1:6)]) {
  
  print(p)
  
  path.output = paste0("4.1_FiguresPopulation/4.1_", p, ".FSTwOutliers.pdf")
  
  png(gsub(".pdf", ".png", path.output), width = 6000, height = 3000, res = 300)
  # pdf(path.output, width = 20, height = 10)
  
  par(mar = c(7.1, 4.1, 1.1, 2.1))
  
  plot(absolute.position, FST.zscores[, p], pch = 20, cex = 0.75, col = col.scaffold,
       xlab = "", ylab = "z-scores (population-specific FST)", xaxt = "n", yaxt = "n",
       ylim = c(floor(min(FST.zscores[, c(7:15)])), ceiling(max(FST.zscores[, c(7:15)]))))
  points(absolute.position[FST.outliers[, p] == 1], FST.zscores[FST.outliers[, p] == 1, p],
         pch = 20, cex = 0.75, col = "#002654")
  
  axis(1, at = thick.mark, labels = gsub("Super-Scaffold_", "", unique(FST.zscores$scaffold)),
       las = 2)
  axis(2, at = seq(-17, 10, 2))

  dev.off()
  
}
