library(tidyverse)
library(readxl)
library(matrixStats)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(42)
dir.create("results", showWarnings = FALSE)

clean_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("MG$", "", x)
  gsub("[^A-Z0-9]", "", x)
}

read_matrix <- function(path) {
  df <- read_excel(path)
  m <- as.matrix(df[, -1])
  rownames(m) <- df[[1]]
  colnames(m) <- clean_id(colnames(m))
  m
}

spearman_ic50 <- function(expr, ic50) {
  keep <- intersect(colnames(expr), names(ic50))
  expr <- expr[, keep]
  ic50 <- ic50[keep]
  s <- apply(expr, 1, function(g) {
    ct <- suppressWarnings(cor.test(g, ic50, method = "spearman", exact = FALSE))
    c(ct$estimate, ct$p.value)
  })
  tibble(gene = rownames(expr), rho = s[1, ], p = s[2, ]) %>%
    mutate(fdr = p.adjust(p, "BH"))
}

volcano <- function(df, title) {
  ggplot(df, aes(rho, -log10(p), colour = fdr <= 0.05)) +
    geom_point(alpha = 0.6, size = 1.6) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#B5613F"), guide = "none") +
    labs(title = title, x = expression(Spearman~rho), y = expression(-log[10]~p)) +
    theme_minimal(base_size = 12)
}

ora_gobp <- function(genes, universe) {
  enrichGO(genes, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = universe,
           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE)
}

save_dotplot <- function(ego, path, title) {
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write_csv(as.data.frame(ego), sub("\\.png$", ".csv", path))
    ggsave(path, dotplot(ego, showCategory = 15) + ggtitle(title),
           width = 8, height = 6, dpi = 300)
  }
}
