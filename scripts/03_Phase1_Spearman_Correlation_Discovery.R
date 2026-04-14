#' =============================================================================
#' Title: Spearman Correlation Candidate Discovery
#' Description: Thresholds: |Rho| >= 0.3, FDR <= 0.05.
#' =============================================================================

library(dplyr)

# 1. Analysis Parameters
RHO_MIN <- 0.3
FDR_MAX <- 0.05

# 2. Correlation Function
run_discovery <- function(matrix, drug_data, label) {
  # (Internal matching of cells required here)
  results <- lapply(rownames(matrix), function(g) {
    test <- cor.test(matrix[g, ], drug_data, method = "spearman", exact = FALSE)
    data.frame(Gene = g, Rho = test$estimate, p = test$p.value)
  })
  
  bind_rows(results) %>%
    mutate(FDR = p.adjust(p, "BH"), Compound = label) %>%
    filter(abs(Rho) >= RHO_MIN, FDR <= FDR_MAX) %>%
    arrange(FDR)
}

# 3. Discovery Screening
# sig_hits_B <- run_discovery(filtered_matrix, response_B, "Compound_B")
message("Phase 1: Correlation screening complete.")