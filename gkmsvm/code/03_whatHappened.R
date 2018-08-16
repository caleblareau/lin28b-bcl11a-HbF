library(BuenColors)
library(dplyr)

importOne <- function() read.table(paste0("../variantPredictions/LIN28B_mendelianRare.out"))

pred_df <- importOne()
