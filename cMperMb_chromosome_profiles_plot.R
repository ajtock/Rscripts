#!/applications/R/R-3.3.2/bin/Rscript

# Plot chromosome-scale cM/Mb profiles using windows
# (stepSize should be 1/10th of window size)

# Usage:
# cMperMb_chromosome_profiles_plot.R 1Mb 10Mb 10-Mb

library(segmentSeq)

args <- commandArgs(trailingOnly = TRUE)
# Define step size and window size
# (stepSize should be 1/10th of window size)
stepName <- as.character(args[1])
winName <- as.character(args[2])
winName2 <- as.character(args[3])

datDir <- "/projects/ajt200/wheat/Wei/cMperMb/"
inDir <- "/projects/ajt200/wheat/Wei/cMperMb/genomeProfiles/"
plotDir <- "/projects/ajt200/wheat/Wei/cMperMb/genomeProfiles/plots/"

dat <- read.table(paste0(datDir, "map.HEI10F2_wRef_opti_1e_8_subset_noNA.txt"),
                  header = T, colClasses = c("NULL", NA, NA, NA))
dat$chr <- paste0("chr", dat$chr)
colnames(dat) <- c("chr", "physical_pos", "genetic_pos")
chrIndex <- unique(dat$chr)
print(chrIndex)
load(file = paste0(inDir, "decimal_suffix_chromosome_lengths.RData"))
print(chrLens)

tabList <- lapply(seq_along(chrIndex), function(x) {
  read.table(file = paste0(inDir, "cMperMb_chromosome_profile_", chrIndex[x],
                           "_win", winName, "_step", stepName, ".txt"))
})
tabList_noNA <- lapply(seq_along(chrIndex), function(x) {
  read.table(file = paste0(inDir, "cMperMb_chromosome_profile_", chrIndex[x],
                           "_win", winName, "_step", stepName, "_noNA.txt"))$x
})

# Function to plot chromosome-scale cM/Mb profiles
cMperMbChrPlot <- function(xplot, dat, dat_noNA, datXlabel) {
  plot(xplot, dat, type = "l", lwd = 1.5, col = "red",
       ylim = c(min(dat_noNA), max(dat_noNA)),
       xlab = "", ylab = "")
  axis(side = 1, lwd.tick = 1.5)
  mtext(side = 1, line = 2.25, cex = 1, text = datXlabel)
  axis(side = 2, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = paste0("cM/Mb (", winName2, " sliding windows)"))
  box(lwd = 1.5)
}

pdf(paste0(plotDir,
           "cMperMb_chromosome_profile_plots_win", winName, "_step", stepName,
           ".pdf"), height = 72, width = 15)
par(mfcol = c(24, 1))
par(mar = c(3.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))
lapply(seq_along(tabList), function(x) {
  cMperMbChrPlot(xplot = tabList[[x]]$windows,
                 dat = tabList[[x]]$cMperMb,
                 dat_noNA = tabList_noNA[[x]],
                 datXlabel = chrIndex[x])
})
dev.off()

