# Cross-platform integration: HGCC vs organoid Spearman rho (Compound B)
source("R/setup.R")

apriori <- c("BRCA1", "BRCA2", "BLM", "RAD51", "FANCD2", "PARP1", "PLK1")

hgcc <- read_csv("results/discovery_CompoundB_spearman.csv", show_col_types = FALSE) %>%
  filter(p < 0.01) %>%
  select(gene, rho_hgcc = rho)

org <- read_csv("results/validation_organoid_spearman.csv", show_col_types = FALSE) %>%
  select(gene, rho_org = rho)

m <- inner_join(hgcc, org, by = "gene") %>%
  mutate(apriori = gene %in% apriori,
         conserved = rho_hgcc < 0 & rho_org < 0)

write_csv(m, "results/integration_cross_platform_rho.csv")
write_csv(filter(m, conserved), "results/conserved_sensitivity_genes.csv")

ggplot(m, aes(rho_hgcc, rho_org)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(aes(colour = apriori), alpha = 0.6, size = 2) +
  scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#B5613F"), guide = "none") +
  geom_text_repel(data = filter(m, apriori), aes(label = gene), size = 3.5) +
  labs(x = expression(HGCC~rho), y = expression(Organoid~rho)) +
  theme_minimal(base_size = 12)

ggsave("results/Fig5b_quadrant.png", width = 6.5, height = 5.5, dpi = 300)
