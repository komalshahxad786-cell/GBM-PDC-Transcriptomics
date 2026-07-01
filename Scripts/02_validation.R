# Validation cohort: 16 patient-derived breast cancer organoids, Compound B
source("R/setup.R")

expr <- read_matrix("data/validation/organoid_expression.xlsx")

ic50 <- read_excel("data/validation/organoid_ic50.xlsx") %>%
  rename_with(~ gsub("IC50 \\(uM\\)", "IC50", .x)) %>%
  mutate(model = clean_id(model)) %>%
  group_by(model) %>%
  summarise(IC50 = median(IC50, na.rm = TRUE), .groups = "drop")

resp <- setNames(log10(ic50$IC50), ic50$model)
cor_org <- spearman_ic50(expr, resp)
write_csv(cor_org, "results/validation_organoid_spearman.csv")

ggsave("results/Fig5a_volcano_organoid.png", volcano(cor_org, "Organoid (Compound B)"),
       width = 6, height = 5, dpi = 300)

ranks <- setNames(-cor_org$rho, cor_org$gene)
ranks <- sort(ranks[!is.na(ranks) & !duplicated(names(ranks))], decreasing = TRUE)

gsea <- gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
              minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.10,
              pAdjustMethod = "BH", verbose = FALSE)

if (!is.null(gsea) && nrow(as.data.frame(gsea)) > 0) {
  write_csv(as.data.frame(gsea), "results/GSEA_GOBP_organoid.csv")
  ggsave("results/Fig6_GSEA_organoid.png",
         dotplot(gsea, showCategory = 15) + ggtitle("Organoid GSEA (GO BP)"),
         width = 8, height = 6, dpi = 300)
}
