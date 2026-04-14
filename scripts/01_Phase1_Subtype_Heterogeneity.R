#' =============================================================================
#' Title: Subtype Heterogeneity Analysis
#' Project: Genomic Determinants of PDC Response in GBM
#' Author: Komal Shahzad
#' Description: Validates phenotypic baseline and evaluates drug response 
#'              distribution across consolidated GBM molecular subtypes.
#' =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)

# 1. Data Loading and Anonymization
setwd("C:/GBM_PDC_Project")
file_xlsx <- "data/raw/Drug_Data_2026.xlsx"

ic50_raw <- read_excel(file_xlsx, sheet = "IC50_Dotmatics")
auc_raw  <- read_excel(file_xlsx, sheet = "AUC_IC50_GraphPad")

# Standardize Cell IDs (Remove 'MG' suffix) and Anonymize Compounds
ic50_clean <- ic50_raw %>%
  rename(Cells = Cells, Compound = Compound, IC50_uM = `IC50 (uM)`) %>%
  mutate(Cells = gsub("MG$", "", as.character(Cells), ignore.case = TRUE),
         Compound = case_when(
           grepl("C-01", Compound) ~ "Compound_A",
           grepl("C-02", Compound) ~ "Compound_B",
           TRUE ~ Compound
         ),
         IC50_uM = as.numeric(IC50_uM)) %>%
  filter(!is.na(IC50_uM))

# 2. Clinical Integration & Subtype Consolidation (PN/NL)
# We group Proneural and Neural as they represent a shared differentiation state
metadata <- auc_raw %>%
  transmute(Cells = gsub("MG$", "", as.character(Cell), ignore.case = TRUE),
            Subtype = as.character(Subtype)) %>%
  distinct()

ic50_final <- ic50_clean %>%
  left_join(metadata, by = "Cells") %>%
  mutate(Subtype_Group = case_when(
    Subtype %in% c("PN", "NL", "Proneural", "Neural") ~ "PN/NL",
    Subtype %in% c("CL", "Classical") ~ "CL",
    Subtype %in% c("MS", "Mesenchymal") ~ "MS",
    TRUE ~ Subtype
  ))

# 3. Figure 1: Subtype Distribution Plot
# Colors: CL=#E41A1C, MS=#377EB8, PN/NL=#4DAF4A
subtype_colors <- c("CL" = "#E41A1C", "MS" = "#377EB8", "PN/NL" = "#4DAF4A")

p1 <- ggplot(ic50_final %>% filter(!is.na(Subtype_Group)), 
             aes(x = Subtype_Group, y = IC50_uM, fill = Subtype_Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1.5) +
  scale_y_log10() +
  scale_fill_manual(values = subtype_colors) +
  facet_wrap(~ Compound) +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_text(face = "bold")) +
  labs(title = "Subtype-Specific Efficacy Profiles",
       x = "GBM Molecular Subtype", y = "IC50 (uM, log10 scale)")

ggsave("results/Fig1_Subtype_Heterogeneity.png", p1, width = 10, height = 6)
message("Phase 1: Subtype analysis complete.")