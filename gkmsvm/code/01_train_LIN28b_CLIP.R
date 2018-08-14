library(gkmSVM)
library(GenomicRanges)
library(diffloop)

if(FALSE){
  gr <- addchr(makeGRangesFromDataFrame(read.table("../../output/LIN28B-IDR-stringent.bed", header = TRUE)))
  write.table(data.frame(gr)[width(gr) < 150,c(1,2,3)], file = "../selected_LIN28B.bed", 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}


do_gkmSVM <- function(){
  
  pop <- "LIN28B"
  genNullSeqs(paste0("../selected_LIN28B.bed"),
              nMaxTrials=10,xfold=1,
              genomeVersion='hg19',
              outputPosFastaFN=paste0('../fasta/',pop,'-positive.fa'),
              outputBedFN=paste0('../fasta/',pop,'-negative.bed'),
              outputNegFastaFN=paste0('../fasta/',pop,'-negative.fa'))
  
  gkmsvm_kernel(paste0('../fasta/',pop,'-positive.fa'),
                paste0('../fasta/',pop,'-negative.fa'),
                paste0('../kernel/', pop, ".kernel.out"))
  
  gkmsvm_trainCV(kernelfn = paste0('../kernel/', pop, ".kernel.out"),
                 posfn = paste0('../fasta/',pop,'-positive.fa'),
                 negfn = paste0('../fasta/',pop,'-negative.fa'), 
                 svmfnprfx=paste0('../kernel/', pop),
                 outputCVpredfn=paste0('../kernel/', pop, ".cvPred.out"),
                 outputROCfn=paste0('../kernel/', pop, ".roc.out"))
  
  gkmsvm_classify('../fasta/nr10mers.fa',
                  svmfnprfx=paste0('../kernel/', pop),
                  paste0('../kernel/', pop, ".weights.10mer.out"))
  pop
  
}

do_gkmSVM()
