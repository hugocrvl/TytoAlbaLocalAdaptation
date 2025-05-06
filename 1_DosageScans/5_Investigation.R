

###############################

########### LOCALLY ###########

###############################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../0_Functions.R", chdir = T)

###############################

# Figure for the IT-PT comparison

###############################

WZA = read.table("../2_WZA/3.6_WZAfinal.txt", header = T)
pairwise.removedNA = read.table("5.2_pairwise.ITPT.DKIS.removedNA.txt",
                                header = T)
coordinates = pairwise.removedNA[, c(1:6)]

middle.position = get.middle.position(coordinates$start, coordinates$end)
absolute.position = get.absolute.position(middle.position, coordinates$scaffold)
thick.mark = get.thick.mark(middle.position, coordinates$scaffold)

pdf("5.3.1_pairwise.ITPT.pdf", width = 20, height = 10)
par(mar = c(7.1, 4.1, 1.1, 2.1))

plot(absolute.position, pairwise.removedNA$IT.PT, xaxt = "n", cex = 0.75,
     col = get.scaffold.colour(coordinates$scaffold), pch = 20, 
     xlab = "", ylab = "Mean kinship between IT and PT individuals")
axis(1, thick.mark, labels = gsub("Super-Scaffold_", "", unique(coordinates$scaffold)),
     las = 2)

dev.off()

###############################

# For the IT-PT pair, comparison
# with the WZA signal

###############################

pdf("5.3.2_pairwise.ITPT.sc22.pdf", width = 15, height = 7.5)

par(oma=c(3,0,0,0))
zones=matrix(c(1,2), ncol=1, byrow=TRUE)
layout(zones, widths=c(1), heights=c(1/2, 1/2))

par(mar = c(0, 5.1, 1.1, 2.1))
plot(absolute.position[coordinates$scaffold == "Super-Scaffold_22"],
     pairwise.removedNA$IT.PT[coordinates$scaffold == "Super-Scaffold_22"],
     pch = 20, col = "#999999", xlab = "", xaxt = "n",
     ylab = "Mean kinship between IT and PT individuals ")

par(mar = c(0, 5.1, 1.1, 2.1))
plot(absolute.position[coordinates$scaffold == "Super-Scaffold_22"],
     -log10(WZA$Z_pVal)[coordinates$scaffold == "Super-Scaffold_22"],
     pch = 20, col = "#999999", xlab = "", ylab = "-log10(p-value)")

dev.off()

pdf("5.3.3_pairwise.ITPT.sc22.correlation.pdf", width = 10, height = 10)

plot(pairwise.removedNA$IT.PT[coordinates$scaffold == "Super-Scaffold_22"],
     -log10(WZA$Z_pVal)[coordinates$scaffold == "Super-Scaffold_22"],
     pch = 20, col = "#999999",
     xlab = "Mean kinship between IT and PT individuals",
     ylab = "-log10(p-value)")

dev.off()

###############################

# Figure for the DK-IS comparison

###############################

pdf("5.4.1_pairwise.DKIS.pdf", width = 20, height = 10)
par(mar = c(7.1,4.1, 1.1, 2.1))

plot(absolute.position, pairwise.removedNA$DK.IS, xaxt = "n", cex = 0.75,
     col = get.scaffold.colour(coordinates$scaffold), pch = 20,
     xlab = "", ylab = "Mean kinship between DK and IS individuals")
axis(1, thick.mark, labels = gsub("Super-Scaffold_", "", unique(coordinates$scaffold)),
     las = 2)

dev.off()

###############################

# For the DK-IS pair, comparison
# with the WZA signal

###############################

pdf("5.4.2_pairwise.DKIS.sc14.pdf", width = 15, height = 7.5)

par(oma=c(3,0,0,0))
zones=matrix(c(1,2), ncol=1, byrow=TRUE)
layout(zones, widths=c(1), heights=c(1/2, 1/2))

par(mar = c(0, 5.1, 1.1, 2.1))
plot(absolute.position[coordinates$scaffold == "Super-Scaffold_14"],
     pairwise.removedNA$DK.IS[coordinates$scaffold == "Super-Scaffold_14"],
     pch = 20, col = "#EAEAEA", xlab = "", xaxt = "n",
     ylab = "Mean kinship between DK and IS individuals")

par(mar = c(0, 5.1, 1.1, 2.1))
plot(absolute.position[coordinates$scaffold == "Super-Scaffold_14"],
     -log10(WZA$Z_pVal)[coordinates$scaffold == "Super-Scaffold_14"],
     pch = 20, col = "#EAEAEA", xlab = "", ylab = "-log10(p-value)")

dev.off()

pdf("5.4.3_pairwise.DKIS.sc14.correlation.pdf", width = 10, height = 10)

plot(pairwise.removedNA$DK.IS[coordinates$scaffold == "Super-Scaffold_14"],
     -log10(WZA$Z_pVal)[coordinates$scaffold == "Super-Scaffold_14"],
     pch = 20, col = "#EAEAEA", 
     xlab = "Mean kinship between DK and IS individuals",
     ylab = "-log10(p-value)")

dev.off()
