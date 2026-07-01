# Transcriptomic Determinants of PDC Sensitivity in Patient-Derived GBM Models

## About

This repository holds the analysis code from my MSc thesis in Molecular Biotechnology
(Systems Biology) at the University of Skövde, carried out at Oncopeptides AB in Stockholm.

The project asks a practical question: can a tumour's gene expression tell us in advance
which tumours will respond to a peptide-drug conjugate (PDC)? Working only from existing
transcriptomic and drug-response data — no new experiments — I correlated gene expression
with drug sensitivity across 48 patient-derived glioblastoma cell cultures (HGCC), took the
signal to the pathway level, and then tested whether it held up in an independent set of 16
patient-derived breast cancer organoids.

The two PDCs converged on the same biology: sensitivity tracked with cell-cycle,
chromosome-segregation and DNA-replication/repair activity rather than with anything
compound-specific, and the pattern reproduced in the organoids. The results are
correlative and hypothesis-generating — candidate determinants of response, not validated
biomarkers.

The two compounds are kept anonymised as **Compound A** and **Compound B**. The code is
here; the drug-response data are proprietary to Oncopeptides and are not included.

I'm a systems biologist with a wet-lab background, interested in turning existing data into
testable biology.
— Komal Shahzad

## Pipeline

**Discovery (HGCC, 2D).** Log2 microarray expression from 48 GBM cell cultures, matched to
IC₅₀ for both compounds. Expression is filtered to the 3,000 most variable genes (MAD),
each gene is correlated with IC₅₀ (Spearman), and p-values are corrected with
Benjamini–Hochberg. Sensitivity- and resistance-associated gene sets (uncorrected P < 0.01)
are passed to GO Biological Process over-representation analysis (ORA).

**Validation (organoids, 3D).** The 150-gene organoid panel is correlated with
log₁₀ IC₅₀ (Compound B), and the list ranked by −ρ is analysed by GSEA against GO BP.

**Integration.** HGCC and organoid Spearman coefficients are compared in a quadrant plot;
seven *a priori* DNA-repair / cell-cycle genes are highlighted, and conserved sensitivity
determinants are defined as ρ < 0 in both cohorts.

IC₅₀ is an inverse measure of sensitivity, so **ρ < 0 marks sensitivity** and ρ > 0 marks
resistance.

## Figures reproduced

| Script | Outputs |
|---|---|
| `01_discovery.R` | Fig 3a/b volcano plots · Fig 4a–c GO BP ORA dotplots · Table 1 (Compound A resistance genes, ρ ≥ 0.3 & FDR ≤ 0.05) |
| `02_validation.R` | Fig 5a organoid volcano · Fig 6 GSEA dotplot |
| `03_integration.R` | Fig 5b cross-platform quadrant · conserved sensitivity gene list |

## Structure

```
R/
├── setup.R           packages and shared functions
├── 01_discovery.R
├── 02_validation.R
├── 03_integration.R
└── run_all.R
data/                 inputs (not tracked — see below)
results/              tables and figures written here
```

## Data

Code only; the input data are not distributed. HGCC expression is public (hgcc.se / GEO),
but the Compound A/B IC₅₀ and the organoid data are proprietary to Oncopeptides AB. Place
your own files in `data/`:

```
data/discovery/hgcc_expression.xlsx     gene symbol + one column per cell line
data/discovery/hgcc_ic50.xlsx           cell_line, Compound_A, Compound_B
data/validation/organoid_expression.xlsx gene symbol + one column per model
data/validation/organoid_ic50.xlsx       model, IC50
```

Adjust the compound / column names in the scripts if yours differ.

## Requirements

R 4.5.1. CRAN: `tidyverse`, `readxl`, `matrixStats`, `ggrepel`. Bioconductor:
`clusterProfiler`, `org.Hs.eg.db`.

```r
install.packages(c("tidyverse", "readxl", "matrixStats", "ggrepel"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
```

## Run

From the repository root:

```r
source("R/run_all.R")
```
