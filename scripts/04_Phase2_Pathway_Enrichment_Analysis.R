#' =============================================================================
#' Title: GO Pathway Enrichment Analysis
#' Description: Maps candidates (|Rho| >= 0.3, FDR <= 0.05) to GO Biological Processes.
#' =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Processing Significant Hits
sig_genes <- sig_hits_B$Gene
entrez_ids <- bitr(sig_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID

# 2. Over-Representation Analysis (ORA)
ego <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", 
                pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

message("Phase 2: Biological pathway enrichment complete.")