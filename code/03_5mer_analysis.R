library(data.table)
library(dplyr)
library(GenomicRanges)
library(BuenColors)
library(diffloop)
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SummarizedExperiment)

"%ni%" <- Negate("%in%")

peaks <- bedToGRanges("../consensus/lin28b_consensus.bed")
counts <- data.matrix(read.table("../consensus/chromVAR_LIN28B.counts.tsv", header = TRUE))

SE <- SummarizedExperiment::SummarizedExperiment(
  rowRanges = peaks, 
  colData = data.frame(samples = c("LIN28B", "PURA")), 
  assays = list(counts = cbind(rowSums(counts[,c(1,2)]), rowSums(counts[,c(3,4)]))))

if(FALSE){
  SE <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = peaks, 
    colData = data.frame(samples = colnames(counts)), 
    assays = list(counts = counts))
}

# Compute deviation scores for 5-mers
SE <- filterPeaks(SE)
kmer_ix <- matchKmers(5, SE, genome = BSgenome.Hsapiens.UCSC.hg19)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)
dev <- computeDeviations(object = SE, annotations = kmer_ix)

df_out <- data.frame(Rep1 = assays(dev)[["z"]][,1],
                     Rep2 = assays(dev)[["z"]][,2], 
                     fiveMer = names(assays(dev)[["z"]][,2]))

df_out %>% arrange(desc(Rep1)) %>% filter(fiveMer == "TCTCC")
