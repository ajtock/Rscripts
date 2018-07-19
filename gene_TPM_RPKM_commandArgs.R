#!/applications/R/R-3.3.2/bin/Rscript

###########################################################################################################
# Calculate gene RPM (reads per million mapped reads), RPKM (reads per kilobase per million mapped reads) #
# and TPM (transcripts per kilobase per million mapped reads) for a given library                         #
###########################################################################################################

# Usage:
# Rscript gene_TPM_RPKM_commandArgs.R ../../WT_RNAseq_meiocyte_Rep1_SRR4204534_STAR_mapped_both_sort.bam WT_RNAseq_meiocyte_Rep1_SRR4204534 SE

source("/projects/ajt200/Rfunctions/RPKMandTPMcalc.R")
library(genomation)
library(ShortRead)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
bamFile <- args[1]
libName <- args[2]
libType <- as.character(args[3])
outDir <- paste0(dirname(bamFile), "/coverage/gene_TPM_RPKM/")

# Load genes as GRanges object
representative_genes_uniq <- system("ls /projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt",
                                    intern=T)
genesGR <- readGeneric(representative_genes_uniq,
                       header=TRUE,
                       strand=4,
                       meta.col=list(gene_model=5))
print("******genes******")
print(genesGR)

# Run RPKMandTPMcalcSE() or RPKMandTPMcalcPE() on the library
myFunc <- paste0("RPKMandTPMcalc", libType)
do.call(myFunc, list(bamFile = bamFile,
                     libName = libName,
                     features = genesGR,
                     featureName = "genes",
                     outDir = outDir))


