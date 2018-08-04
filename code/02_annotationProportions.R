library(data.table)
library(dplyr)
library(GenomicRanges)
library(BuenColors)
library(diffloop)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

"%ni%" <- Negate("%in%")

# Import annotation
x01 <- bedToGRanges("../annotations/UCSC_3primeUTR.bed")
x02 <- bedToGRanges("../annotations/UCSC_Exons.bed")
x03 <- bedToGRanges("../annotations/UCSC_5primeUTR.bed")
x04 <- bedToGRanges("../annotations/UCSC_Introns.bed")

# Take a sample and 
gr_t <- addchr(makeGRangesFromDataFrame(read.table("../output/LIN28B-Rep1.bed", header = TRUE)))
seq_t <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr_t)
motif <- as.numeric((vcountPattern("GGAGA", seq_t) + vcountPattern("TCTCC", seq_t)) > 0)

# Do overlaps
ov_1 <- findOverlaps(gr_t, x01)
ov_2 <- findOverlaps(gr_t, x02)
ov_3 <- findOverlaps(gr_t, x03)
ov_4 <- findOverlaps(gr_t, x04)

# Classify each variant
class <- ifelse(1:length(gr_t) %in% queryHits(ov_1), "3UTR",
                ifelse(1:length(gr_t) %in% queryHits(ov_2), "Exon",
                       ifelse(1:length(gr_t) %in% queryHits(ov_3), "5UTR",
                              ifelse(1:length(gr_t) %in% queryHits(ov_4), "Intron", "other"))))

data.frame(class, motif) %>%
  group_by(motif, class) %>% 
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
