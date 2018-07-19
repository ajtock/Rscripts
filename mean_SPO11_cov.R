#!/applications/R/R-3.3.2/bin/Rscript

##################################################################################
# Make mean normalised coverage tab-separated bed file for 3 SPO11-1-oligos reps #
##################################################################################

covDir <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/"
SPO11_rep1 <- read.table(file = paste0(covDir, "WT_SPO11-oligo_RPI1_norm_allchrs_coverage_coord_tab.bed"))
SPO11_rep2 <- read.table(file = paste0(covDir, "WT_SPO11-oligo_RPI3_norm_allchrs_coverage_coord_tab.bed"))
SPO11_rep3 <- read.table(file = paste0(covDir, "WT_SPO11-oligo_RPI8_norm_allchrs_coverage_coord_tab.bed"))

SPO11_allreps <- cbind(SPO11_rep1, SPO11_rep2[,4], SPO11_rep3[,4])
SPO11_allreps_mean <- cbind(SPO11_allreps,
  (SPO11_allreps[,4] + SPO11_allreps[,5] + SPO11_allreps[,6]) / 3)
SPO11_mean <- cbind(SPO11_allreps_mean[,1:3], SPO11_allreps_mean[,7])
write.table(SPO11_mean, file = paste0(covDir, "WT_SPO11-oligo_meanAllReps_norm_allchrs_coverage_coord_tab.bed"),
  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

