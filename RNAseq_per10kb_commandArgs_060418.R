#!/applications/R/R-3.3.2/bin/Rscript

#############################################################################################
# Calculate coverage values in 10-kb windows for generating chromosome-scale plots          #
#############################################################################################

# Usage on hydrogen node7
# csmit -m 20G -c 1 "Rscript RNAseq_per10kb_commandArgs_060418.R ../WT_RNAseq_Rep1_norm_allchrs_coverage_coord_tab.bed WT_RNAseq_Rep1"

library(segmentSeq)

args <- commandArgs(trailingOnly = TRUE)
covFile <- args[1]
ChIPname <- args[2]
outDir <- paste0(dirname(covFile), "/genomeProfiles/") 

covDat <- read.table(covFile)

# genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# pericentromeric regions are as defined in Supplemental Table S26 of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

windows <- c(10000)
winNames <- c("10kb")

#################################
# make cumulative genomes       #
#################################

sumchr <- cumsum(c(0, chrLens))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

# calculate mean coverage values within windows and log2 transform
for(s in 1:length(windows)) {
  print(s)
  covWins <- windows[s]
  print(winNames[s])
  cumChIPwinDat <- NULL
    for(i in 1:5) {
      print(i)
      chrChIPinput <- covDat[covDat[,1] == chrs[i],]
      ChIPcov <- chrChIPinput[,4]
      covCoords <- seq(1, length(ChIPcov), by = 1)
      covIRcoords <- IRanges(start = covCoords, width = 1)
      covGRcoords <- GRanges(seqnames = chrs[i], strand = "+", ranges = covIRcoords)
      seqWindows <- seq(1, chrLens[i], by = covWins)
      seqWindows <- c(seqWindows, chrLens[i])
      cumWindows <- seqWindows + sumchr[i]
      windowsIRanges <- IRanges(start = seqWindows, width = covWins)
      windowsGRanges <- GRanges(seqnames = chrs[i], strand = "+", ranges = windowsIRanges)
      overlaps <- getOverlaps(windowsGRanges, covGRcoords, whichOverlaps = TRUE)
      ChIPcovWinVals <- sapply(overlaps, function(x) mean(ChIPcov[x])+1) #+1 is offset to avoid infinite values
      log2ChIPcovWinVals <- log2(ChIPcovWinVals)
      ChIPwinDat <- cbind(cumWindows, log2ChIPcovWinVals)
      cumChIPwinDat <- rbind(cumChIPwinDat, ChIPwinDat)
      write.table(ChIPwinDat, file = paste0(outDir, "log2_", ChIPname, "_chr", i, "_norm_coverage_", winNames[s], ".txt"))
    }
    write.table(cumChIPwinDat, file = paste0(outDir, "log2_", ChIPname, "_genome_norm_coverage_", winNames[s], ".txt"))
}

chrlog2trans <- lapply(1:5, function(i) {
                  read.table(file = paste0(outDir, "log2_", ChIPname, "_chr", i, "_norm_coverage_", winNames[1], ".txt"))
                })

# smooth coverage values with MA filter
test <- seq(1, 1000, by = 1)
j = 100
ma <- rep(1, test[j])/test[j]

filt_log2trans <- NULL
filt_log2trans_noNA <- NULL  
  for(i in 1:5) {
    filt_chrlog2trans <- stats::filter(chrlog2trans[i][[1]][,2], ma)
    which_na <- which(is.na(filt_chrlog2trans) == TRUE)
    left_na <- which_na[which(which_na < 100)]
    left_val <- filt_chrlog2trans[left_na[length(left_na)]+1]
    filt_chrlog2trans[left_na] <- left_val
    right_na <- which_na[which(which_na > 100)]
    right_val <- filt_chrlog2trans[right_na[1]-1]
    filt_chrlog2trans[right_na] <- right_val
    filt_chrlog2trans_noNA <- filt_chrlog2trans[!is.na(filt_chrlog2trans)]
    filt_chrlog2trans <- cbind(chrlog2trans[i][[1]][,1], filt_chrlog2trans)
    write.table(filt_chrlog2trans, file = paste0(outDir, "filt_log2_", ChIPname, "_chr", i, "_norm_coverage_", winNames[1], ".txt"))
    write.table(filt_chrlog2trans_noNA, file = paste0(outDir, "filt_noNA_log2_", ChIPname, "_chr", i, "_norm_coverage_", winNames[1], ".txt"))
    filt_log2trans <- rbind(filt_log2trans, filt_chrlog2trans)
    filt_log2trans_noNA <- c(filt_log2trans_noNA, filt_chrlog2trans_noNA)
  }
write.table(filt_log2trans, file = paste0(outDir, "filt_log2_", ChIPname, "_genome_norm_coverage_", winNames[1], ".txt"))
write.table(filt_log2trans_noNA, file = paste0(outDir, "filt_noNA_log2_", ChIPname, "_genome_norm_coverage_", winNames[1], ".txt"))

