#' =============================================================================
#' Title: Mechanistic Discovery Plots
#' Description: Volcano plots and targeted barplots. Highlight threshold: p < 0.01.
#' =============================================================================

library(ggplot2)
library(stringr)

# 1. Targeted Sensitivity Barplot (Compound_A)
# Highlight: Firebrick for Cell Cycle/Replication, Gray for others
df_plot <- ego_results %>%
  mutate(Highlight = grepl("cell cycle|division|replication", Description, ignore.case = TRUE))

p5 <- ggplot(df_plot, aes(x = Count, y = reorder(Description, Count), fill = Highlight)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("gray70", "firebrick")) +
  theme_minimal() +
  labs(title = "Compound_A: Proliferation Prerequisite Drivers",
       subtitle = "Selection Criteria: p < 0.01", x = "Gene Count", y = "")

ggsave("results/Fig5_Sensitivity_Barplot.png", p5, width = 8, height = 6)