# GBM-PDC-Transcriptomics
Bioinformatic pipeline investigating transcriptomic drivers of peptide-drug conjugate (PDC) sensitivity and resistance in Glioblastoma (GBM).
# Genomic Determinants of PDC Response in Glioblastoma

**Author:** Komal Shahzad  
**Contact:** komal.shahxad786@gmail.com

## 🔬 Project Overview & Objectives
This project integrates human glioblastoma cell culture (HGCC) gene expression data with continuous drug response profiles (AUC/IC50) to identify biomarkers of sensitivity. The pipeline specifically focuses on a comparative analysis of two novel Peptide-Drug Conjugates (PDCs): Compound A and Compound B(aliases used to protect proprietary pharmaceutical confidentiality). By mapping transcriptomic signatures, the ultimate goal is to define the precise Mechanism of Action (MoA) and identify targetable clinical vulnerabilities.

## 🛠️ Computational Workflow
* **Data Harmonization:** Matched patient-derived cell models across HGCC expression matrices and proprietary drug datasets.
* **Phase 1 (Broad Screen):** Analyzed a broader panel of compounds to establish a global phenotypic baseline and identify macro-level sensitivity patterns.
* **Phase 2 (Mechanism Deep Dive):** Focused specifically on Compound A and and Compound B, utilizing continuous Spearman rank correlations and GO/Reactome pathway enrichment to isolate the shared mechanism of action.
* **Phase 3 (Validation):** Cross-referenced 2D cell line signatures against 3D patient-derived organoid models (Crown Biosciences).

## 🧬 Key Findings: The "Two-Hit" Mechanism of Action

### 1. The "Master Switch": Proliferation vs. Differentiation
The efficacy of these compounds is strongly influenced by the baseline transcriptional state of the tumor cells. 
* **Sensitivity Drivers:** Maximum sensitivity is strictly driven by high proliferation. Pathway enrichment (GO/Reactome) confirms that genes regulating the cell cycle, mitosis, and DNA replication consistently predict compound efficacy.
* **Resistance Drivers:** Resistant cells exhibit highly differentiated phenotypes, relying heavily on extracellular matrix (ECM) organization, cell migration, and mesenchymal transitions to escape drug toxicity.

### 2. Conserved Genomic Signatures
Comparative hierarchical clustering revealed that Compound A and Compound B exhibit identical genomic response profiles. This striking symmetry—where genes driving sensitivity for Compound A do the exact same for Compound B—confirms that both compounds functionally target the same biological axis: **Replication-Dependent DNA Damage.**

### 3. The "Two-Hit" Model of Sensitivity
The data supports a definitive "Two-Hit" model dictating PDC efficacy:
* **Hit 1 (The Prerequisite):** The cell must be actively dividing. This baseline proliferation creates the necessary replication stress required for the drug to induce double-strand DNA breaks.
* **Hit 2 (The Determinant):** The cell must be deficient in Homologous Recombination (HR) DNA repair. Cells with low repair capacities cannot withstand the induced stress and undergo apoptosis.

## 🧫 Phase 3 Findings: 3D Organoid Validation (Crown Biosciences)
To ensure these findings were not mere 2D cell culture artifacts, the transcriptomic signature was validated against complex 3D patient-derived organoids treated with Compound B (acting as the representative compound for the PDC class).

* **Cross-Model Conservation:** A targeted quadrant analysis comparing HGCC (2D) versus CrownBio (3D) correlation scores revealed highly conserved biology. Genes driving sensitivity and resistance in flat cell lines maintained their distinct mechanistic roles in 3D organoids.
* **Mechanism Confirmed:** The 3D organoid Volcano and Enrichment plots perfectly mirrored the Phase 2 findings, successfully confirming that **DNA Repair Deficiency and Replication Stress** are robust, translatable biomarkers for PDC efficacy in complex biological systems.

## 📂 Repository Structure
* `scripts/`: Contains the modular R pipeline for variance filtration, continuous modeling, and pathway enrichment.
* `results/`: High-resolution outputs of Volcano plots, Cross-Model Quadrant plots, Heatmaps, and GO/Reactome Enrichment visualizations.
* `data/`: *Raw transcriptomic matrices, CrownBio organoid data, and proprietary IC50/AUC viability data are strictly withheld to comply with pharmaceutical confidentiality agreements. (See `.gitignore`)*
