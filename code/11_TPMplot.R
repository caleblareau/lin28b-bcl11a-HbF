library(data.table)
library(dplyr)
library(BuenColors)

ad1 <- data.frame(fread("../data/Ad1-kallisto/abundance.tsv"))
ad2 <- data.frame(fread("../data/Ad2-kallisto/abundance.tsv"))
cb1 <- data.frame(fread("../data/CB1-kallisto/abundance.tsv"))
cb2 <- data.frame(fread("../data/CB2-kallisto/abundance.tsv"))

df <- data.frame(
  transcript = ab1$target_id,
  AD1 = ad1$tpm,
  AD2 = ad2$tpm,
  CB1 = cb1$tpm,
  CB2 = cb2$tpm
)

BCL11A <- c("ENST00000492272.5",
            "ENST00000489516.6", 
            "ENST00000489183.1",
            "ENST00000477659.1",
            "ENST00000409351.4",
            "ENST00000359629.9",
            "ENST00000358510.5",
            "ENST00000356842.8",
            "ENST00000335712.10")

df %>% filter(transcript %in% BCL11A) %>%
  reshape2::melt(id.vars = "transcript") %>%
  mutate(condition = ifelse(substring(variable, 1, 1) == "A", "Adult", "Newborn")) %>%
  group_by(transcript, condition) %>% summarize(mean = mean(log2(value + 1)), sd = sd(log2(value + 1))) -> summary_df

p1 <- ggplot(summary_df, aes(x = transcript, y = mean, fill = condition)) +
  geom_bar(position=position_dodge(), stat="identity", color = "black", width = 0.4) +
  scale_fill_manual(values = c("lightgrey", "gray30")) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd),
                width=.2, position=position_dodge(.4)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "", y = paste0("log2 TPM")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.position = "bottom")
ggsave(p1, file = "../output/TPM-Kallisto.pdf", width = 4, height = 3)


df %>% 
  mutate(meanNB = (CB1 + CB2)/2, 
         meanAD = (AD1 + AD2)/2) -> allDF

