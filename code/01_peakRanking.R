library(data.table)
library(dplyr)
library(GenomicRanges)
library(BuenColors)
library(Biostrings)
library(diffloop)
library(BSgenome.Hsapiens.UCSC.hg19)

"%ni%" <- Negate("%in%")

LIN28B_binding <- makeGRangesFromDataFrame(data.frame(
  chr = "2", start = 60688540, end = 60688629
))

importPeaks <- function(file){
  df <- data.frame(fread(file))
  x <- unique(makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE))
  x <- x[seqnames(x) %in% c(as.character(1:22), "X")]
  return(x)
}

lin_rna <- importPeaks("../peaks/LIN28B_noDup_noLambda_peaks.narrowPeak")
lin_rrna <- importPeaks("../peaks/LIN28B_rRNA_noDup_noLambda_peaks.narrowPeak")
p1 <- importPeaks("../peaks/PURA-Abcam_noDup_noLambda_peaks.narrowPeak")
p2 <- importPeaks("../peaks/PURA-Bethyl_noDup_noLambda_peaks.narrowPeak")

lin28b_RNA_NoDupNoLambda_go <- lin_rna[1:length(lin_rna) %ni%  queryHits(findOverlaps(lin_rna, c(p1, p2)))]
seq_lin28b_RNA_NoDupNoLambda_go <- getSeq(BSgenome.Hsapiens.UCSC.hg19, addchr(lin28b_RNA_NoDupNoLambda_go))

motif_RNA <- as.numeric((vcountPattern("GGAGA", seq_lin28b_RNA_NoDupNoLambda_go) + vcountPattern("TCTCC", seq_lin28b_RNA_NoDupNoLambda_go)) > 0)
sumdf_RNA <- data.frame(lin28b_RNA_NoDupNoLambda_go, motif = motif_RNA, BCL11A = 1:length(lin28b_RNA_NoDupNoLambda_go) %in%
                      queryHits(findOverlaps(lin28b_RNA_NoDupNoLambda_go, LIN28B_binding)))

sumdf_RNA %>% arrange(desc(V9)) %>% filter(motif == 1) %>% 
  mutate(percentile_log10Q = 1:n()/n(), rank = 1:n(), total = n()) %>% filter(BCL11A) -> rankDF_RNA


lin28b_rRNA_NoDupNoLambda_go <- lin_rrna[1:length(lin_rrna) %ni%  queryHits(findOverlaps(lin_rrna, c(p1, p2)))]
seq_lin28b_rRNA_NoDupNoLambda_go <- getSeq(BSgenome.Hsapiens.UCSC.hg19, addchr(lin28b_rRNA_NoDupNoLambda_go))
motif_rRNA <- as.numeric((vcountPattern("GGAGA", seq_lin28b_rRNA_NoDupNoLambda_go) + vcountPattern("TCTCC", seq_lin28b_rRNA_NoDupNoLambda_go)) > 0)

sumdf_rRNA <- data.frame(lin28b_rRNA_NoDupNoLambda_go, motif = motif_rRNA, BCL11A = 1:length(lin28b_rRNA_NoDupNoLambda_go) %in%
                          queryHits(findOverlaps(lin28b_rRNA_NoDupNoLambda_go, LIN28B_binding)))

sumdf_rRNA %>% arrange(desc(V9))  %>% 
  mutate(percentile_log10Q = 1:n()/n(), rank = 1:n(), total = n()) %>% filter(BCL11A) -> rankDF_rRNA

rankDF_rRNA
rankDF_RNA

sumdf_rRNA %>% arrange(desc(V9)) %>% mutate(rank = 1:n()) %>%
  arrange(BCL11A) %>% 
  ggplot(aes(x = rank, y = V9, color = BCL11A)) + geom_point(size = 0.2) +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "Rank Sorted Peaks - LIN28B Rep1", y = "-log10 Q-value") +
  scale_color_manual(values = c("grey", "black")) + theme(legend.position = "none") -> P1

sumdf_RNA %>% arrange(desc(V9)) %>% mutate(rank = 1:n()) %>%
  arrange(BCL11A) %>% 
  ggplot(aes(x = rank, y = V9, color = BCL11A)) + geom_point(size = 0.2) +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "Rank Sorted Peaks - LIN28B Rep2", y = "-log10 Q-value") +
  scale_color_manual(values = c("grey", "black")) + theme(legend.position = "none") -> P2

cowplot::ggsave(cowplot::plot_grid(P1, P2, nrow = 1), file = "../output/plots/rankSortedPeaks.pdf",  width = 4, height = 2)

out_rep1 <- sumdf_rRNA[,c(1,2,3,9,10,11)]; colnames(out_rep1) <- c("chr", "start", "end", "fold-change", "log10P", "log10Q")
write.table(out_rep1, file = "../output/LIN28B-Rep1.bed", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

out_rep2 <- sumdf_RNA[,c(1,2,3,9,10,11)]; colnames(out_rep2) <- c("chr", "start", "end", "fold-change", "log10P", "log10Q")
write.table(out_rep2, file = "../output/LIN28B-Rep2.bed", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
            

   

