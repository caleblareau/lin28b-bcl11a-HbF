library(seqinr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gkmSVM)


# Function to compute delta 
do_gkmsvm_delta <- function(){
  population <- "LIN28B"
  gkmsvm_delta('../fasta/human-FromMouseMismatch.fasta','../fasta/mouse-FromMouseMismatch.fasta',
               svmfnprfx=paste0('../kernel/', population),
               paste0('../variantPredictions/',population,'_mendelianRare.out'))
}

do_gkmsvm_delta()



