#!/applications/R/R-3.3.2/bin/Rscript

# Generate chromosome-scale cM/Mb profiles using windows
# (stepSize should be 1/10th of window size)

# Usage:
# cMperMb_chromosome_profiles.R 100000 100kb 1Mb

library(segmentSeq)

args <- commandArgs(trailingOnly = TRUE)
# Define step size and window size
# (stepSize should be 1/10th of window size)
stepSize <- as.numeric(args[1])
stepName <- as.character(args[2])
winName <- as.character(args[3])

inDir <- "/projects/ajt200/wheat/Wei/cMperMb/"
outDir <- "/projects/ajt200/wheat/Wei/cMperMb/genomeProfiles/"

dat <- read.table(paste0(inDir, "map.HEI10F2_wRef_opti_1e_8_subset_noNA.txt"),
                  header = T, colClasses = c("NULL", NA, NA, NA))
dat$chr <- paste0("chr", dat$chr)
colnames(dat) <- c("chr", "physical_pos", "genetic_pos")

chrIndex <- unique(dat$chr)

physical_delta <- NULL
genetic_delta <- NULL
for(i in 1:length(chrIndex)) {
  datChr <- dat[dat$chr == chrIndex[i],]
  physical_delta_chr <- as.numeric(sapply(seq(from = 0, to = dim(datChr)[1]-1),
    function(x) {
      datChr$physical_pos[x+1]-datChr$physical_pos[x]
    }))
  physical_delta <- c(physical_delta, physical_delta_chr)
  genetic_delta_chr <- as.numeric(sapply(seq(from = 0, to = dim(datChr)[1]-1),
    function(x) {
      datChr$genetic_pos[x+1]-datChr$genetic_pos[x]
    }))
  genetic_delta <- c(genetic_delta, genetic_delta_chr)
}

dat <- cbind(dat, physical_delta, genetic_delta,
             genetic_delta/(physical_delta/1000000))
colnames(dat) <- c("chr", "physical_pos", "genetic_pos",
                   "physical_delta", "genetic_delta", "cMperMb")

# For each chromosome with decimal suffix,
# define chromosome length as max SNP coordinate
chrIndex <- unique(dat$chr)
chrLens <- NULL
for(i in 1:length(chrIndex)) {
  datChr <- dat[dat$chr == chrIndex[i],]
  chrLens <- c(chrLens, max(datChr$physical_pos, na.rm = T))
}

# For each chromosome, define chromosome length as max SNP coordinate
allProbes <- read.table(paste0(inDir, "allprobegenome.blast"))
chrIndex2 <- levels(allProbes$V2)[-length(levels(allProbes$V2))]
chrLens2 <- NULL
for(i in 1:length(chrIndex2)) {
  allProbesChr <- allProbes[allProbes$V2 == chrIndex2[i],]
  chrLens2 <- c(chrLens2, max(allProbesChr$V9))
}

print(chrIndex[1:5])
#[1] "chr1A.1" "chr1B.1" "chr1D.1" "chr2A.1" "chr2B.4"
print(chrIndex2[1:5])
#[1] "chr1A" "chr1B" "chr1D" "chr2A" "chr2B"
chrLens[1:5] <- chrLens2[1:5]

print(chrIndex[7:12])
#[1] "chr2D.1" "chr3A.1" "chr3B.1" "chr3D.2" "chr4A"   "chr4B.1"
print(chrIndex2[6:11])
#[1] "chr2D" "chr3A" "chr3B" "chr3D" "chr4A" "chr4B"
chrLens[7:12] <- chrLens2[6:11]

print(chrIndex[13])
#[1] "chr5A.1"
print(chrIndex2[13])
#[1] "chr5A"
chrLens[13] <- chrLens2[13]

print(chrIndex[15:16])
#[1] "chr5B.1" "chr5D.1"
print(chrIndex2[14:15])
#[1] "chr5B" "chr5D"
chrLens[15:16] <- chrLens2[14:15]

print(chrIndex[18:21])
#[1] "chr6A"   "chr6B"   "chr6D.1" "chr7A.1"
print(chrIndex2[16:19])
#[1] "chr6A" "chr6B" "chr6D" "chr7A"
chrLens[18:21] <- chrLens2[16:19]

print(chrIndex[23:24])
#[1] "chr7B.1" "chr7D.3"
print(chrIndex2[20:21])
#[1] "chr7B" "chr7D"
chrLens[23:24] <- chrLens2[20:21]

print("Chromosome names")
print(chrIndex)
print("Chromosome lengths")
print(chrLens)

save(chrLens, file = paste0(outDir, "decimal_suffix_chromosome_lengths.RData"))

library(doParallel)
registerDoParallel(cores = length(chrIndex))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(i = 1:length(chrIndex)) %dopar% {
  # Define steps
  stepStarts <- seq(from = 1, to = chrLens[i], by = stepSize)
  #stepEnds <- c(stepStarts[2:length(stepStarts)]-1, chrLens[i])
  #stepGR <- GRanges(seqnames = chrIndex[i],
  #                  ranges = IRanges(start = stepStarts, end = stepEnds))
  winGR <- GRanges()
  for(j in 1:(length(stepStarts)-10)) {
    winGRtmp <- GRanges(seqnames = chrIndex[i],
                        ranges = IRanges(start = stepStarts[j],
                                         end = stepStarts[j+10]-1))
    winGR <- append(winGR, winGRtmp)
  }
  winGR <- append(winGR,
                  GRanges(seqnames = chrIndex[i],
                          ranges = IRanges(start = stepStarts[length(stepStarts)],
                                           end = chrLens[i])))
  # Define SNP intervals
  datChr <- dat[dat$chr == chrIndex[i],]
  snpGR <- GRanges()
  for(k in 1:(dim(datChr)[1]-1)) {
    snpGRtmp <- GRanges(seqnames = chrIndex[i],
                        ranges = IRanges(start = datChr$physical_pos[k],
                                         end = datChr$physical_pos[k+1]),
                        cMperMb = datChr$cMperMb[k+1])
    snpGR <- append(snpGR, snpGRtmp)
  }
  overlaps <- getOverlaps(winGR, snpGR, whichOverlaps = T)
  winVals <- sapply(overlaps,
                    function(x) mean(snpGR$cMperMb[x]))
  winDat <- cbind(round(start(winGR)+((end(winGR)-start(winGR))/2)), winVals)
  colnames(winDat) <- c("windows", "cMperMb")
  
  # Replace NA values with leftmost or rightmost value
  which_na <- which(is.na(winDat[,2]) == TRUE)
  left_na <- which_na[which(which_na < dim(winDat)[1]/2)]
  left_val <- winDat[,2][left_na[length(left_na)]+1]
  winDat[,2][left_na] <- left_val
  right_na <- which_na[which(which_na > dim(winDat)[1]/2)]
  right_val <- winDat[,2][right_na[1]-1]
  winDat[,2][right_na] <- right_val
  winDat_noNA <- winDat[,2][!is.na(winDat[,2])] 
  write.table(winDat,
              file = paste0(outDir, "cMperMb_chromosome_profile_", chrIndex[i],
                            "_win", winName, "_step", stepName, ".txt"))
  write.table(winDat_noNA,
              file = paste0(outDir, "cMperMb_chromosome_profile_", chrIndex[i],
                            "_win", winName, "_step", stepName, "_noNA.txt"))
}


