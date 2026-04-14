# Genomic Determinants of PDC Response in Glioblastoma
**Author:** Komal Shahzad  
**Contact:** komal.shahxad786@gmail.com

## 🔬 Project Overview 
This repository contains a complete, modular R pipeline designed to investigate the transcriptomic drivers of peptide-drug conjugate (PDC) sensitivity and resistance in Glioblastoma (GBM). 

By integrating human glioblastoma cell culture (HGCC) RNA-seq matrices with continuous drug response profiles (IC50/AUC), this pipeline identifies actionable biomarkers and elucidates the mechanism of action (MoA) for two novel compounds (aliases: Compound A and Compound B). Findings derived from 2D HGCC screens were subsequently validated against 3D patient-derived organoid models (Crown Biosciences).

## 🛠️ Technical Stack & Dependencies
* **Language:** R (v4.x)
* **Data Wrangling:** `dplyr`, `tidyr`, `purrr`, `readxl`, `stringr`
* **Pre-processing:** `matrixStats` (Median Absolute Deviation filtration)
* **Statistical Analysis:** Base R `cor.test` (Spearman Rank, BH-adjusted FDR)
* **Enrichment & Annotation:** `clusterProfiler`, `org.Hs.eg.db`
* **Visualization:** `ggplot2`, `ggrepel`
* ## 📊 Methodology & Statistical Parameters
To ensure strict reproducibility and robust biomarker selection, the following statistical thresholds and methodologies were standardized across the pipeline:

* **Variance Filtration:** Pre-processed HGCC RNA-seq matrices were log2-transformed and filtered using Median Absolute Deviation (MAD). The pipeline isolated the top **3,000 highly variable genes (HVGs)** to reduce computational noise and focus on transcriptionally active targets.
* **Continuous Drug Modeling (Phase 1):** Biomarker discovery relied on continuous Spearman rank correlations between HVG expression and drug response metrics. Candidate genes were strictly thresholded at an effect size of **|Rho| >= 0.3** and a significance of **FDR <= 0.05** (Benjamini-Hochberg correction).
* **Pathway Enrichment:** Over-Representation Analysis (ORA) mapped significant signatures to GO Biological Processes, utilizing a **BH-adjusted p-value cutoff of 0.05**.
* **3D Organoid Validation (Phase 2):** To account for the inherent variance in complex patient-derived 3D models, cross-model correlation validation utilized a targeted significance threshold of **|Rho| >= 0.3** and an unadjusted **p-value < 0.01**.

## ⚙️ Computational Pipeline & Usage
The analysis is structured sequentially. Scripts are designed to be run in numeric order. *(Note: Raw expression matrices and proprietary compound data are withheld per confidentiality agreements).*

* **`01_Phase1_Subtype_Heterogeneity.R`**: Ingests raw IC50/AUC data, standardizes cell IDs, maps models to consolidated GBM molecular subtypes (PN/NL, CL, MS), and visualizes baseline phenotypic distribution.
* **`02_Phase1_Variance_Filtration.R`**: Pre-processes HGCC RNA-seq data via log2 transformation and Median Absolute Deviation (MAD) filtration to isolate the top 3,000 highly variable genes (HVGs) for downstream modeling.
* **`03_Phase1_Spearman_Correlation_Discovery.R`**: Executes continuous transcriptomic modeling. Computes Spearman rank correlations between HVG expression and drug response, applying strict thresholds (`|Rho| >= 0.3`, `FDR <= 0.05` via Benjamini-Hochberg).
* **`04_Phase2_Pathway_Enrichment_Analysis.R`**: Performs Over-Representation Analysis (ORA) on significant gene signatures against Gene Ontology (GO) Biological Processes.
* **`05_Phase2_Mechanistic_Visualizations.R`**: Generates targeted statistical visualizations (Volcano plots, barplots) highlighting pathway enrichment (e.g., cell cycle, DNA replication drivers).
* **`07_Phase2_3D_Organoid_Validation.R`**: A robust validation pipeline that ingests 3D CrownBio organoid data, harmonizes passage metrics, and executes cross-model correlation mapping (2D HGCC vs. 3D Organoid). Outputs Quadrant validation plots for targeted DNA repair panels.

## 🧬 Key Biological Findings

### 1. The "Two-Hit" Mechanism of Action
The transcriptomic pipeline revealed a highly conserved "Two-Hit" prerequisite for PDC efficacy:
* **Hit 1 (Proliferation):** Sensitivity is strictly driven by high baseline proliferation, creating replication stress. GO enrichment confirmed cell cycle, mitosis, and DNA replication as primary efficacy predictors.
* **Hit 2 (DNA Repair Deficiency):** Cells lacking robust Homologous Recombination (HR) and DNA damage repair capacity undergo apoptosis upon PDC exposure. Conversely, resistant models exhibit differentiated phenotypes relying on ECM organization and mesenchymal transitions.

### 2. 3D Organoid Translation (CrownBio)
To ensure the identified mechanisms were not 2D culture artifacts, the 2D HGCC genomic signatures were mapped against complex 3D patient-derived organoids treated with Compound B.
* **Cross-Model Conservation:** Quadrant analysis of Spearman rho values demonstrated strict conservation of sensitivity drivers.
* **Target Validation:** Critical DNA repair and replication genes—including **BRCA1, BRCA2, BLM, RAD51, PARP1, and PLK1**—maintained their distinct mechanistic roles in the 3D organoid environment, validating them as robust, translatable biomarkers for this PDC class.

## 📂 Repository Structure

├── results/          
│   ├── CompoundA_Sensitivity_Enrichment.xlsx
│   ├── CompoundB_Sensitivity_Enrichment.xlsx
│   ├── Selected_Variable_Genes.xlsx
│   ├── Fig1_Subtype_Efficacy_Distribution.png
│   ├── Fig2_CompoundA_Discovery_Volcano.jpg
│   ├── Fig3_CompoundB_Discovery_Volcano.jpg
│   ├── Fig4_Top_Hits_Correlation_Heatmap.png
│   ├── Fig5_CompoundA_Proliferation_Barplot.png
│   ├── Fig5b_CompoundB_DNA Repair_Barplot.png
│   ├── Fig7a_Organoid_Validation_Volcano.png
│   ├── Fig7b_2D_vs_3D_Quadrant_Validation.png
│   ├── Fig8a_CompoundB_Sensitivity_GO_Dotplot_organoidValidation.png
│   └── Fig8b_CompoundB_Resistance_GO_Dotplot_organoidValidation.png
├── scripts/
│   ├── 01_Phase1_Subtype_Heterogeneity.R
│   ├── 02_Phase1_Variance_Filtration.R
│   ├── 03_Phase1_Spearman_Correlation_Discovery.R
│   ├── 04_Phase2_Pathway_Enrichment_Analysis.R
│   ├── 05_Phase2_Mechanistic_Visualizations.R
│   └── 07_Phase2_3D_Organoid_Validation.R
├── .gitignore
├── LICENSE
└── README.md
