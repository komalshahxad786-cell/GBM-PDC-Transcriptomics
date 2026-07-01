# Discovery cohort: 48 HGCC GBM cell cultures, Compound A and Compound B
source("R/setup.R")

expr <- read_matrix("data/discovery/hgcc_expression.xlsx")

ic50 <- read_excel("data/discovery/hgcc_ic50.xlsx") %>%
  rename_with(~ gsub("IC50 \\(uM\\)", "IC50", .x)) %>%
  mutate(cell_line = clean_id(cell_line)) %>%
  group_by(cell_line) %>%
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop")

expr <- expr[, colnames(expr) %in% ic50$cell_line]
expr <- expr[order(rowMads(expr), decreasing = TRUE)[1:3000], ]

resp <- function(cmp) setNames(ic50[[cmp]], ic50$cell_line)
cor_a <- spearman_ic50(expr, resp("Compound_A"))
cor_b <- spearman_ic50(expr, resp("Compound_B"))
write_csv(cor_a, "results/discovery_CompoundA_spearman.csv")
write_csv(cor_b, "results/discovery_CompoundB_spearman.csv")

ggsave("results/Fig3a_volcano_CompoundA.png", volcano(cor_a, "Compound A"),
       width = 6, height = 5, dpi = 300)
ggsave("results/Fig3b_volcano_CompoundB.png", volcano(cor_b, "Compound B"),
       width = 6, height = 5, dpi = 300)

universe <- rownames(expr)
save_dotplot(ora_gobp(filter(cor_a, p < 0.01, rho < 0)$gene, universe),
             "results/Fig4a_ORA_CompoundA_sensitivity.png", "Compound A sensitivity")
save_dotplot(ora_gobp(filter(cor_b, p < 0.01, rho < 0)$gene, universe),
             "results/Fig4b_ORA_CompoundB_sensitivity.png", "Compound B sensitivity")
save_dotplot(ora_gobp(filter(cor_b, p < 0.01, rho > 0)$gene, universe),
             "results/Fig4c_ORA_CompoundB_resistance.png", "Compound B resistance")

cor_a %>%
  filter(rho >= 0.3, fdr <= 0.05) %>%
  arrange(fdr) %>%
  write_csv("results/Table1_CompoundA_resistance_genes.csv")
