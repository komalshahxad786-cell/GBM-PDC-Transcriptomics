#' =============================================================================
#' Title: 3D Organoid Validation Pipeline
#' Project: GBM PDC Transcriptomic Analysis
#' Author: Komal Shahzad
#' Description: Validates transcriptomic mechanisms in 3D patient-derived models.
#'              Outputs Volcano, Cross-Model Quadrant, and GO Enrichment plots.
#' Thresholds: |Rho| >= 0.3, p < 0.01.
#' =============================================================================

# 0) Environment Setup
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)

setwd("C:/GBM_PDC_Project")
data_dir <- "data/raw"
results_dir <- "results"
if(!dir.exists(results_dir)) dir.create(results_dir)

# 1) Data Ingestion
message("Loading CrownBio 3D Organoid datasets...")
drug_raw <- read_excel(file.path(data_dir, "Crown_organoid_IC50_all.xlsx"))
expr_raw <- read_excel(file.path(data_dir, "Gene_Expression_RNAseq_Crown.xlsx"))

# 2) Data Harmonization (The "Passage Eraser")
clean_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

drug_final <- drug_raw %>%
  rename(Model_ID = `Organoid Name`, Compound = `Test Article`, IC50 = `Relative IC50 (µM)`) %>%
  filter(grepl("Compound_B", Compound, ignore.case = TRUE)) %>%
  mutate(
    Model_ID = clean_id(Model_ID),
    IC50 = as.numeric(gsub(",", ".", IC50)),
    Compound = "Compound_B" 
  ) %>%
  drop_na(IC50)

expr_final <- expr_raw %>%
  rename(Model_ID = Model, Gene = Gene, Value = Value) %>%
  mutate(
    Model_ID = clean_id(Model_ID),
    Model_ID = gsub("P[0-9]+$", "", Model_ID), # Strip Passage suffixes
    Value = as.numeric(Value)
  ) %>%
  group_by(Model_ID, Gene) %>%
  summarize(Value = mean(Value, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Gene, values_from = Value)

merged_data <- inner_join(drug_final, expr_final, by = "Model_ID")
message("Success: Merged ", nrow(merged_data), " models for validation.")

# 3) Full-Transcriptome Spearman Correlation
message("Executing Spearman correlation analysis...")
gene_list <- setdiff(names(merged_data), c("Model_ID", "Compound", "IC50"))

cor_results <- map_dfr(gene_list, function(g) {
  vec <- as.numeric(merged_data[[g]])
  if(sum(!is.na(vec)) < 3 || var(vec, na.rm=TRUE) == 0) return(NULL)
  
  res <- cor.test(vec, log10(merged_data$IC50), method = "spearman", exact = FALSE)
  data.frame(Gene = g, Rho = as.numeric(res$estimate), P_value = res$p.value)
}) %>%
  mutate(NegLogP = -log10(P_value)) %>%
  arrange(P_value)

# ==============================================================================
# VISUALIZATION 1: ORGANOID VOLCANO PLOT
# ==============================================================================
p_volcano <- ggplot(cor_results, aes(x = Rho, y = NegLogP)) +
  geom_point(color = "grey85", alpha = 0.4, size = 1.5) +
  geom_point(data = filter(cor_results, P_value < 0.01 & Rho <= -0.3), 
             color = "#2b8cbe", alpha = 0.7, size = 2) +
  geom_point(data = filter(cor_results, P_value < 0.01 & Rho >= 0.3), 
             color = "#e34a33", alpha = 0.7, size = 2) +
  geom_text_repel(data = head(cor_results, 15), aes(label = Gene), 
                  size = 3.5, fontface = "bold", box.padding = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted", color = "grey50") +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot: Crown Organoids Compound B",
       subtitle = paste0("N = ", nrow(merged_data), " Models | |Rho| >= 0.3, p < 0.01"),
       x = "Spearman Rho (Negative = Sensitivity)", y = "-log10(p-value)")

ggsave(file.path(results_dir, "Fig6_CompoundB_Organoid_Validation_Volcano.png"), p_volcano, width = 8, height = 6)

# ==============================================================================
# VISUALIZATION 2: TARGETED QUADRANT PLOT (2D vs 3D)
# ==============================================================================
if(exists("sig_hits_B")) {
  quadrant_data <- inner_join(
    sig_hits_B %>% select(Gene, HGCC_Rho = Rho),
    cor_results %>% select(Gene, Organoid_Rho = Rho),
    by = "Gene"
  )
  
  target_genes <- c("BRCA1", "BRCA2", "BLM", "RAD51", "FANCD2", "PARP1", "PLK1")
  
  p_quadrant <- ggplot(quadrant_data, aes(x = HGCC_Rho, y = Organoid_Rho)) +
    geom_point(color = "grey60", alpha = 0.5, size = 2) +
    geom_point(data = filter(quadrant_data, Gene %in% target_genes), 
               color = "#CC0000", size = 4) +
    geom_text_repel(data = filter(quadrant_data, Gene %in% target_genes), 
                    aes(label = Gene), color = "#CC0000", fontface = "bold", size = 5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    annotate("text", x = -0.3, y = -0.4, label = "Conserved Sensitivity", 
             color = "#5C76E6", size = 4.5, hjust = 0) +
    theme_minimal(base_size = 14) +
    labs(title = "Compound B Sensitivity: HGCC vs Crown Organoids",
         subtitle = "Targeted Validation of DNA Repair/Cell Cycle Panel",
         x = "HGCC Correlation (Spearman rho)",
         y = "Crown Organoids Correlation (Spearman rho)")
  
  ggsave(file.path(results_dir, "Fig7_CompoundB_Quadrant_Validation.png"), p_quadrant, width = 9, height = 7)
} else {
  message("Skipping Quadrant Plot: Phase 1 HGCC data not found in environment.")
}

# ==============================================================================
# VISUALIZATION 3: GO ENRICHMENT (SENSITIVITY & RESISTANCE)
# ==============================================================================
message("Performing GO Pathway Enrichment for Organoid Hits...")

# Isolate strict hits
sens_genes <- cor_results %>% filter(Rho <= -0.3, P_value < 0.01) %>% pull(Gene)
res_genes  <- cor_results %>% filter(Rho >= 0.3, P_value < 0.01) %>% pull(Gene)

# Helper function to run enrichment and plot
run_go_enrichment <- function(gene_list, plot_title, filename) {
  if(length(gene_list) < 5) {
    message("Not enough genes to run GO for: ", plot_title)
    return(NULL)
  }
  
  entrez <- tryCatch(
    bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
    error = function(e) return(NULL)
  )
  
  if(is.null(entrez)) return(NULL)
  
  ego <- enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP", 
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  
  if(!is.null(ego) && nrow(ego) > 0) {
    p_dot <- dotplot(ego, showCategory = 15) + 
      ggtitle(plot_title) + 
      theme_bw()
    
    ggsave(file.path(results_dir, filename), p_dot, width = 9, height = 7)
    return(p_dot)
  }
  return(NULL)
}

# Generate Fig 8a: Sensitivity GO
run_go_enrichment(
  gene_list = sens_genes, 
  plot_title = "Compound B Organoid Sensitivity (GO BP)", 
  filename = "Fig8a_CompoundB_Organoid_Sensitivity_GO.png"
)

# Generate Fig 8b: Resistance GO
run_go_enrichment(
  gene_list = res_genes, 
  plot_title = "Compound B Organoid Resistance (GO BP)", 
  filename = "Fig8b_CompoundB_Organoid_Resistance_GO.png"
)

write.xlsx(cor_results, file.path(results_dir, "Organoid_Validation_Full_Results.xlsx"))
message("Phase 2: Organoid validation pipeline complete.")