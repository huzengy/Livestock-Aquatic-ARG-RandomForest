# Random Forest analysis of ARG subtypes

R scripts for binary classification of Chongming Island surface-water and adjacent Yangtze River background samples based on ARG subtype abundances (copies/cell).

## Scripts

- `01_RF_feature_selection.R` — screens ARG subtypes using Wilcoxon rank-sum tests, Benjamini–Hochberg FDR adjustment, directional mean abundance, fold change (≥1.5), and prevalence; trains a 1,000-tree Random Forest with `set.seed(123)`; exports the top 20 features ranked by Mean Decrease Gini.
- `02_RF_model_validation.R` — independently repeats the same subtype screening, fits the full-data model, calculates Gini and permutation importance, checks Spearman collinearity, and performs **five repetitions of stratified five-fold cross-validation**. Within each training fold, it repeats candidate screening and top-20 selection before testing on held-out samples. Outputs CV metrics (accuracy, balanced accuracy, AUC), out-of-fold predictions, importance rankings, selection frequencies, and diagnostic tables.

**No plotting or figure-export code is included.**

## Input

An Excel workbook named `ARG_subtype_abundance.xlsx` in the working directory, with:

- `Water`: Chongming Island surface-water samples; first two columns are metadata.
- `Sheet1`: reference samples; the first column contains a group named `Yangtze River`, and the first two columns are metadata.
- Remaining shared, numeric columns: ARG subtype abundances in copies/cell. Missing values are interpreted as zero, matching the original analysis.

The scripts expect 30 island and 23 background-river samples for the validation step. **The input workbook is not included in this repository.** The scripts therefore cannot be run from this repository alone without the abundance matrix in the indicated format.

## Requirements and usage

R packages: `readxl`, `dplyr`, `randomForest` (the validation script only requires `readxl` and `randomForest`).

```r
source("01_RF_feature_selection.R")
source("02_RF_model_validation.R")
```

Results are written to `results/` and `results/validation/`, respectively. Both scripts retain the original model settings, filtering thresholds, and random seeds. The validation script also reports permutation-importance ranks, using the `randomForest` importance implementation. It performs cross-validation independently rather than reading previously calculated metrics from CSV files.

## Scope

This repository contains the authors' Random Forest feature-selection and model-validation code only. Standard metagenomic processing and other statistical analyses are described in the manuscript's Materials and Methods.
