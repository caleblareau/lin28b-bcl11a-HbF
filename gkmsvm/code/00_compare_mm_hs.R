library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(GenomicRanges)
library(seqinr)


strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

# Known coordinates from liftover
mouse <- makeGRangesFromDataFrame(
  data.frame(chr = "chr11", start = 24164037, end = 24164175)
)
human <- makeGRangesFromDataFrame(
  data.frame(chr = "chr2", start = 60688530, end = 60688668)
)

# Sequence specification
mm_seq <- as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, mouse))
hs_seq <- strReverse(chartr("ATGC","TACG",as.character(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, human)))))

lapply(10:129, function(i){
  if(substring(hs_seq, i, i) != substring(mm_seq, i, i)){
    print(i)
    substring(mm_seq, (i-9), (i + 9))
  } else {
    ""
  }
}) -> mm_gkmer

lapply(10:129, function(i){
  if(substring(hs_seq, i, i) != substring(mm_seq, i, i)){
    print(i)
    substring(hs_seq, (i-9), (i + 9))
  } else {
    ""
  }
}) -> hs_gkmer

names <- paste0("mut", as.character(1:6))

if(FALSE){
  write.fasta(hs_gkmer[c(17, 38, 54, 62, 74, 104)-9], names, "../fasta/human-FromMouseMismatch.fasta",
              open = "w", nbchar = 60)
  write.fasta(mm_gkmer[c(17, 38, 54, 62, 74, 104)-9], names, "../fasta/mouse-FromMouseMismatch.fasta",
              open = "w", nbchar = 60)
}