library(data.table)
library(dplyr)

df <- fread("GSE107218_CBPB-hg19-counts.txt.gz")
df2 <- df[,c(7:54)]
colnames(df2) <- gsub("-sorted.bam", "", gsub("hisat", "", colnames(df2))) 
mat <- df2 %>% data.frame() %>% data.matrix()
tpm <- round(t(t(mat)/colSums(mat))*1000000, 1)
tpm_out <- data.frame(tpm, gene = df[,1])
write.table(tpm_out, file = "GSE107218_TPM_converted.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)