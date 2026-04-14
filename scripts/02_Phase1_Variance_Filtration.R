#' =============================================================================
#' Title: Expression Matrix Pre-processing & Filtration
#' Description: Normalizes HGCC RNA-seq and filters for highly variable genes.
#' =============================================================================

library(readxl)
library(matrixStats)
library(dplyr)

# 1. Loading HGCC Transcriptome
expr_raw <- read_excel("data/raw/HGCC_expression_all_genes.xlsx")

# 2. Cleaning and Numeric Coercion
expr_mat <- expr_raw %>%
  distinct(Gene.Symbol, .keep_all = TRUE) %>%
  filter(!is.na(Gene.Symbol)) %>%
  select(-Gene.Symbol, -matches("ProbeSetID")) %>%
  as.matrix() %>%
  apply(2, as.numeric)

rownames(expr_mat) <- expr_raw$Gene.Symbol[1:nrow(expr_mat)]
expr_mat[is.na(expr_mat)] <- 0

# 3. Normalization and MAD Filtration
log_expr <- log2(expr_mat + 1)
mad_values <- rowMads(log_expr)

# Retaining Top 3,000 genes for primary correlation screening
filtered_matrix <- log_expr[order(mad_values, decreasing = TRUE)[1:3000], ]

saveRDS(filtered_matrix, "data/processed/Filtered_Expression_Matrix.rds")
message("Phase 1: Variance filtration complete. 3000 HVGs retained.")