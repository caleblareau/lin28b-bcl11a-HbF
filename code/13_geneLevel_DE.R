library(sleuth)
library(BuenColors)
library(dplyr)

base_dir <- "../rna-seqdata/kallisto"
sample_id <- dir(file.path(base_dir))
kal_dirs <- paste0(base_dir, "/", sample_id)

s2c <- data.frame(
 sample = c("Ad1", "Ad2", "CB1", "CB2"),
 condition = c("Ad", "Ad", "CB", "CB"), 
 path = kal_dirs,
 stringsAsFactors = FALSE
)

# Do biomart things to collapse transcripts to gene level counts
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'grch37.ensembl.org')

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g,
                  aggregation_column = 'ens_gene', gene_mode = TRUE, transformation_function = function(x) log2(x + 0.5))

so <- sleuth_fit(so)
so <- sleuth_wt(so, "conditionCB")
results_table2 <- sleuth_results(so, 'conditionCB', test_type = 'wt')
plot_volcano(so, 'conditionCB', test_type = 'wt')

p1 <- ggplot(results_table2, aes(x = b, y = -log10(qval))) +
  geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() + theme(legend.position = "none") +
  scale_y_continuous(expand = c(0.01,0.01))
cowplot::ggsave(p1, file = "../output/plots/volcano_sleuth.pdf", width  = 2.5, height = 2.5)

       