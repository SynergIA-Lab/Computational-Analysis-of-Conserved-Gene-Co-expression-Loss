# Computational Analysis of Conserved Gene Co-expression Loss Reveals Prognostic Regulatory Network Disruption Across Solid Tumors

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)

## Overview

This repository contains the complete analytical pipeline used in the study:

> Ríos-Cadenas M., Segura-Carmona I., López-Fernández A., Gómez-Vela F.A. *Computational Analysis of Conserved Gene Co-expression Loss Reveals Prognostic Regulatory Network Disruption Across Solid Tumors*. Preprint submitted to Computational Biology and Chemistry, 2026.

The pipeline identifies gene co-expression interactions that are strongly preserved in healthy tissue but consistently lost during malignant transformation across four solid tumor types from The Cancer Genome Atlas (TCGA): breast invasive carcinoma (BRCA), lung adenocarcinoma (LUAD), head and neck squamous cell carcinoma (HNSC), and stomach adenocarcinoma (STAD).

---

## Repository Structure

```
├── 01_download_and_stratify.R         # Data acquisition and cohort stratification
├── 02_preprocess_quality_control.R    # Quality control and outlier removal
├── 03_differential_expression_analysis.R  # DEG identification (DESeq2 + edgeR)
├── 04_coexpression_network_inference.R    # Co-expression network construction
├── 05_lost_interactions_consensus.R       # Cross-cohort consensus lost interactions
├── 06_network_evaluation.R                # GeneMANIA-based network validation
├── 07_survival_analysis.R                 # Prognostic validation (Kaplan-Meier, Cox)
└── README.md
```

> **Note on working directories**: Each script contains a `BASE_DIR` or `setwd()` parameter at the top that must be updated to match your local directory structure before execution.

---

## Data Availability

All input data are publicly available from The Cancer Genome Atlas (TCGA) through the [GDC Data Portal](https://portal.gdc.cancer.gov/). The project identifiers used are:

- `TCGA-BRCA` — Breast Invasive Carcinoma
- `TCGA-LUAD` — Lung Adenocarcinoma
- `TCGA-HNSC` — Head and Neck Squamous Cell Carcinoma
- `TCGA-STAD` — Stomach Adenocarcinoma

RNA-seq count data (HTSeq-Counts / STAR-Counts workflow) and clinical metadata were downloaded using the GDC Data Transfer Tool.

The GeneMANIA co-expression reference network used for validation is available at [https://genemania.org/](https://genemania.org/).

---

## Requirements

### R version
R >= 4.0 is required.

### R packages

The following packages are required. Install them via Bioconductor and CRAN:

```r
# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "limma",
  "apeglm"
))

# CRAN packages
install.packages(c(
  "WGCNA",
  "survival",
  "survminer",
  "dplyr",
  "readr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "RColorBrewer",
  "data.table",
  "matrixStats",
  "stringr",
  "tibble",
  "cluster",
  "vegan",
  "gridExtra"
))
```

---

## Pipeline Description

The analytical pipeline consists of 7 sequential steps. Each script must be executed in order as the output of each step serves as the input for the next.

---

### Step 1 — Data Acquisition and Cohort Stratification
**Script**: `01_download_and_stratify.R`

Downloads raw RNA-seq count data (STAR-Counts workflow) from TCGA using the `TCGAbiolinks` package and organizes samples into biologically meaningful subgroups.

**Key operations**:
- Downloads RNA-seq data for the specified TCGA project via the GDC API
- Performs sample deduplication (one sample per patient, prioritized by library size)
- Maps Ensembl gene IDs to gene symbols and collapses duplicates by summation
- Stratifies samples by molecular or anatomical subtype:
  - BRCA → PAM50 subtypes (Basal, Her2, LumA, LumB)
  - LUAD → Anatomical location (Upper lobe, Lower lobe)
  - HNSC → Anatomical subsite (Larynx, Oral Cavity, Oropharynx)
  - STAD → Molecular subtype (CIN, MSI)

**Input**: TCGA project identifier (e.g., `"TCGA-BRCA"`)

**Output** (per subtype directory):
```
counts_renombrado.csv     # Full count matrix (all samples)
counts_primary.csv        # Tumor samples only
counts_control.csv        # Normal tissue samples only
mapeo_muestras.csv        # Sample ID mapping file
```

**Configuration to update**:
```r
projects <- c("TCGA-BRCA")   # Change to desired TCGA project
setwd("/your/working/directory")
```

---

### Step 2 — Preprocessing and Quality Control
**Script**: `02_preprocess_quality_control.R`

Performs comprehensive multi-metric quality control to identify and remove aberrant samples prior to network analysis.

**Key operations**:
- Filters low-abundance genes using `edgeR::filterByExpr` (CPM-based)
- Applies Variance Stabilizing Transformation (VST) via `DESeq2`
- Generates quality control visualizations:
  - Library size barplots
  - PCA plots (VST-normalized)
  - Hierarchical clustering heatmaps (Euclidean, Pearson, and Spearman distances)
  - MDS plots (Pearson and Spearman distances)
- Detects outlier samples through multi-metric consensus (PCA ±3 SD, hierarchical clustering, MDS)
- Validates group separation using PERMANOVA before and after outlier removal
- Exports the final list of samples to remove

**Input**:
```
counts_primary.csv        # From Step 1
counts_control.csv        # From Step 1
```

**Output**:
```
initial_qc/               # QC figures (PDF)
FINAL_samples_to_remove.txt   # List of outlier samples
permanova_result_*.txt    # PERMANOVA results
```

**Configuration to update**:
```r
setwd("/your/path/TCGA-BRCA")
path_csv <- "./by_Subtype/LumB/counts_primary.csv"
ctrl_csv <- "./by_Subtype/LumB/counts_control.csv"
```

---

### Step 3 — Differential Expression Analysis
**Script**: `03_differential_expression_analysis.R`

Identifies differentially expressed genes (DEGs) between tumor and control samples using DESeq2, applying stringent statistical thresholds.

**Key operations**:
- Loads count matrices and removes outlier samples identified in Step 2
- Filters low-expression genes using `edgeR::filterByExpr`
- Runs `DESeq2` differential expression analysis (Tumor vs. Control)
- Applies VST normalization for downstream network analysis
- Identifies DEGs at two fold change thresholds: |log2FC| > 2 and |log2FC| > 3 (FDR < 0.05, Benjamini-Hochberg correction)
- Generates volcano plots and heatmaps of top DEGs

**Input**:
```
counts_primary.csv            # From Step 1
counts_control.csv            # From Step 1
FINAL_samples_to_remove.txt   # From Step 2
```

**Output** (per FC threshold, in `DEGs_analysis/FC2/` and `DEGs_analysis/FC3/`):
```
normal_DEGs_FC*.csv           # Control expression matrix (DEGs only)
tumor_DEGs_FC*.csv            # Tumor expression matrix (DEGs only)
DESeq2_results_complete_FC*.csv   # Full DESeq2 results table
activated_genes_FC*.txt       # Upregulated gene list
repressed_genes_FC*.txt       # Downregulated gene list
all_DEGs_FC*.txt              # Complete DEG list
volcanoPlot_FC*.pdf/.png      # Volcano plots
heatmap_top50_DEGs_FC*.pdf    # Heatmap of top 50 DEGs
```

**Configuration to update**:
```r
setwd("/your/path/TCGA-BRCA")
samples_to_remove_file <- "./by_Subtype/Basal/FINAL_samples_to_remove.txt"
```

---

### Step 4 — Co-expression Network Inference
**Script**: `04_coexpression_network_inference.R`

Constructs separate Pearson correlation-based co-expression networks for tumor and control conditions, restricted to DEGs, and identifies lost interactions using a double-threshold filtering strategy.

**Key operations**:
- Loads DEG expression matrices (tumor and control) from Step 3
- Applies `WGCNA::goodSamplesGenes` for final data quality check
- Computes Pearson correlation matrices independently for tumor and control conditions
- Identifies lost interactions satisfying all three criteria simultaneously:
  - Control correlation ≥ 0.7 (strong co-expression in healthy tissue)
  - Tumor correlation ≤ 0.3 (loss of co-expression in tumor)
  - ΔCorrelation ≥ 0.3 (minimum magnitude of change)
- Exports correlation matrices and lost interaction catalogue

**Input**:
```
normal_DEGs_FC*.csv       # From Step 3
tumor_DEGs_FC*.csv        # From Step 3
```

**Output** (in `resultsWGCNA/FC*/`):
```
lost_interactions.csv         # Full catalogue of lost interactions
correlation_matrix_control.rds    # Control correlation matrix
correlation_matrix_tumor.rds      # Tumor correlation matrix
genes_used.txt                # List of genes included in analysis
analysis_summary.txt          # Summary report
```

**Configuration to update**:
```r
BASE_DIR <- "/your/path/TCGA-BRCA/by_Subtype/LumB"
INPUT_DIR <- file.path(BASE_DIR, "DEGs_analysis/FC3")
```

---

### Step 5 — Cross-cohort Consensus Lost Interactions
**Script**: `05_lost_interactions_consensus.R`

Identifies lost interactions conserved across multiple cancer types through a two-level hierarchical comparison: within-cancer consensus and cross-cancer pairwise intersection.

**Key operations**:
- Loads lost interaction catalogues from all 11 biological contexts (4 BRCA subtypes, 2 LUAD locations, 3 HNSC subsites, 2 STAD subtypes)
- **Project-level consensus**: retains only interactions present in ALL subdatasets within each cancer type (robust to intratumor heterogeneity)
- **Pairwise intersection**: identifies interactions shared between all pairs of cancer types (BRCA-LUAD, BRCA-HNSC, BRCA-STAD, LUAD-HNSC, LUAD-STAD, HNSC-STAD)
- Computes mean ΔCorrelation across datasets for each consensus interaction
- Processes both FC2 and FC3 thresholds

**Input**:
```
lost_interactions.csv     # From Step 4, for each of the 11 subdatasets
```

**Output** (in `_Common_Interactions/FC*/`):
```
01_PROYECTO_TCGA-*.csv    # Project-level consensus interactions per cancer type
02_PAR_TCGA-*_vs_TCGA-*.csv   # Pairwise cross-cancer interactions
*_Resumen.txt             # Summary reports for each comparison
```

**Configuration to update**:
```r
base_path <- "/your/working/directory"
datasets <- list(
  "TCGA-LUAD/by_Tissue/Upper_lobe__lung",
  ...   # Update paths to match your directory structure
)
```

---

### Step 6 — Network Validation Against GeneMANIA
**Script**: `06_network_evaluation.R`

Evaluates the biological validity of identified lost interactions by comparing them against the GeneMANIA gold-standard co-expression database using a confusion matrix framework.

**Key operations**:
- Loads lost interaction catalogues from Step 4
- Loads the GeneMANIA co-expression reference network
- Filters GeneMANIA to the gene universe of each analysis
- Constructs a confusion matrix (TP, FP, FN, TN) over all possible gene pair interactions
- Computes precision, recall, F-score, accuracy, and specificity
- Generates detailed evaluation reports for FC2 and FC3 thresholds
- Produces a comparative summary table across threshold configurations

**Input**:
```
lost_interactions.csv         # From Step 4
red_coexpresion_final.txt     # GeneMANIA co-expression network (downloaded separately)
```

**Output** (in `resultsWGCNA/FC*/`):
```
evaluation_metrics_FC*.csv    # Performance metrics per FC threshold
confusion_matrix_FC*.csv      # Confusion matrix per FC threshold
network_evaluation_report_FC*.txt   # Detailed evaluation report
comparison_summary.csv        # Comparative table across thresholds
```

**Configuration to update**:
```r
CONFIG$base_dir <- "/your/path/TCGA-STAD/by_Molecular_Subtype/MSI"
CONFIG$genemania_file <- "/your/path/red_coexpresion_final.txt"
```

---

### Step 7 — Prognostic Validation and Survival Analysis
**Script**: `07_survival_analysis.R`

Assesses the clinical relevance of conserved lost interactions by linking residual co-expression patterns to patient overall survival using Kaplan-Meier estimation and Cox proportional hazards regression.

**Key operations**:
- Loads cross-cancer pairwise interaction sets from Step 5 (LUAD-BRCA, LUAD-STAD, HNSC-STAD)
- Computes a pairwise co-expression score for each patient: S = z_g1 · z_g2
- Stratifies patients into tertiles (low, medium, high co-expression groups)
- Performs Kaplan-Meier survival analysis with log-rank test
- Performs Cox proportional hazards regression reporting HR with 95% CI for medium and high groups vs. low reference
- Generates Kaplan-Meier curves with risk tables for all significant interactions (log-rank p < 0.05)
- Exports consolidated results table ordered by significance

**Input**:
```
02_PAR_TCGA-*_vs_TCGA-*.csv   # From Step 5 (pairwise interactions)
counts_renombrado.csv         # From Step 1 (expression data per subdataset)
TCGA-*_metadata.tsv           # Clinical metadata from TCGA GDC portal
```

**Output** (in `Survival_Analysis_Results/`):
```
Figures/KM_*.pdf              # Kaplan-Meier curves per interaction and cancer type
Tables/survival_results_all.csv   # Complete survival analysis results
Reports/survival_analysis_summary.txt  # Summary report
```

**Configuration to update**:
```r
CONFIG$base_dir <- "/your/working/directory"
CONFIG$common_interactions_dir <- "_Common_Interactions/your_date/FC2"
```

---

## Execution Order Summary

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 1 | `01_download_and_stratify.R` | TCGA GDC Portal | Count matrices per subtype |
| 2 | `02_preprocess_quality_control.R` | Count matrices | QC figures, outlier list |
| 3 | `03_differential_expression_analysis.R` | Count matrices + outlier list | DEG lists + expression matrices |
| 4 | `04_coexpression_network_inference.R` | DEG expression matrices | Lost interactions per subdataset |
| 5 | `05_lost_interactions_consensus.R` | Lost interactions (all subdatasets) | Consensus lost interactions |
| 6 | `06_network_evaluation.R` | Lost interactions + GeneMANIA | Performance metrics |
| 7 | `07_survival_analysis.R` | Consensus interactions + expression + clinical | KM curves + Cox HR |

> Steps 2–4 must be executed independently for each of the 11 biological contexts (subdatasets). Steps 5–7 are run once using the combined outputs from all subdatasets.

---

## Citation

If you use this code in your research, please cite:

```
Ríos-Cadenas M., Segura-Carmona I., López-Fernández A., Gómez-Vela F.A.
Computational Analysis of Conserved Gene Co-expression Loss Reveals Prognostic 
Regulatory Network Disruption Across Solid Tumors.
Preprint submitted to Computational Biology and Chemistry, 2026.
```

---

## Contact

For questions regarding the code or methodology, please contact:

**Francisco Gómez-Vela, Ph.D.**  
SynergIA Laboratory (SIALAB)  
Universidad Pablo de Olavide, Seville, Spain  
📧 fgomez@upo.es

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
