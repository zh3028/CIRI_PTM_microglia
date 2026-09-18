# PTM Enzyme Co-expression Code in Microglia after CIRI

Analysis code for the manuscript:

> **Three PTM Enzyme Co-expression Modules Define Microglial Functional Heterogeneity after Cerebral Ischemia-Reperfusion Injury**  
> Ying Li, Huimin Li, Meng Zhang*  
> *Corresponding author: Meng Zhang (zh3028@zuaa.zju.edu.cn)*

---

## Overview

This repository contains the R scripts and processed data for single-cell transcriptomic analysis of microglia after cerebral ischemia-reperfusion injury (CIRI). We identified three post-translational modification (PTM) enzyme co-expression modules:

- **M1**: Metabolic stress-associated
- **M2**: Pro-inflammatory-associated
- **M3**: Reparative-associated

These modules define distinct microglial functional states with dynamic spatiotemporal transitions.

---

## Repository Structure

```
.
├── README.md
├── LICENSE
├── .gitignore
├── sessionInfo.txt
├── code/                          # R scripts used for all analyses
└── results/                       # Processed data and analysis outputs
    ├── D1_H_matrix_correct.csv
    ├── D1_W_matrix.csv
    ├── D1_W_matrix_correct.csv
    ├── D1_nmf_full.rds
    ├── D1_nmf_full_correct.rds
    ├── D1_nmf_input_matrix.csv
    ├── D1_nmf_input_matrix_dense.csv
    ├── D2_H_matrix.csv
    ├── D2_H_matrix_loose.csv
    ├── nrun_stability_report.csv
    ├── rank_selection_metrics.csv
    └── tables/
        ├── D1_D5_correlation.csv
        ├── D1_module_1_GO.csv
        ├── D1_module_2_GO.csv
        ├── D1_module_3_GO.csv
        ├── D1_module_markers.csv
        ├── D1_top_genes.txt
        ├── D1_top_markers.csv
        ├── D2_module_proportions.csv
        ├── D3_validation_results.csv
        ├── D4_sex_module_table.csv
        └── ptm_genes_final.txt
```

---

## Data Sources

All raw data are publicly available from the Gene Expression Omnibus (GEO):

| Dataset | GEO Accession | Model | Timepoints | Use |
| :--- | :--- | :--- | :--- | :--- |
| D1 | GSE174574 | tMCAO reperfusion | 24h | Module identification |
| D2 | GSE227651 | tMCAO reperfusion | 1d/3d/7d | Spatiotemporal dynamics |
| D3 | GSE245386 | tMCAO reperfusion | 24h | Independent validation |
| D4 | GSE267240 | Permanent ischemia | — | Sex difference analysis |
| D5 | GSE319237 | Permanent ischemia | 14d | Purity validation (bulk) |

---

## Requirements

- **R** version 4.5.2 (primary) and 4.3.3 (for NMF stability testing)
- Key R packages:
  - `Seurat` (v5.0.3)
  - `NMF` (v0.28)
  - `clusterProfiler` (v4.10.0)
  - `org.Mm.eg.db` (v3.22.0)
  - `Monocle3` (v1.3.1)
  - `ggplot2`, `dplyr`, `tidyr`, `patchwork`

Install dependencies:

```r
install.packages(c("Seurat", "NMF", "clusterProfiler", "org.Mm.eg.db", 
                   "Monocle3", "ggplot2", "dplyr", "tidyr", "patchwork"))
```

---

## Usage

Run the R scripts in the `code/` folder in order. Intermediate and final outputs are saved to the `results/` folder.

```bash
Rscript code/01_preprocessing.R
Rscript code/02_microglia_filtering.R
Rscript code/03_NMF_modules.R
Rscript code/04_GO_enrichment.R
Rscript code/05_spatiotemporal.R
Rscript code/06_validation.R
Rscript code/07_purity_validation.R
```

---

## Key Results

| Figure | Description | Output file |
| :--- | :--- | :--- |
| Figure 1 | Study design flowchart | — |
| Figure 2 | Microglial purity validation | results/tables/D1_D5_correlation.csv |
| Figure 3 | NMF module identification | results/D1_W_matrix.csv |
| Figure 4 | GO enrichment of modules | results/tables/D1_module_*_GO.csv |
| Figure 5 | Spatiotemporal dynamics | results/tables/D2_module_proportions.csv |
| Figure 6 | Independent validation | results/tables/D3_validation_results.csv |
| Figure S1 | Quality control violin plots | results/ |
| Figure S2 | Sex-specific module proportions | results/tables/D4_sex_module_table.csv |

---

## Citation

If you use this code or data in your research, please cite:

> Li Y, Li H, Zhang M. Three PTM Enzyme Co-expression Modules Define Microglial Functional Heterogeneity after Cerebral Ischemia-Reperfusion Injury. *Frontiers in Immunology* (under review).

A preprint is available on bioRxiv:

> Li Y, Li H, Zhang M. A PTM Regulatory Enzyme Co-expression Code Defines Microglial Functional Heterogeneity in Cerebral Ischemia-Reperfusion Injury. *bioRxiv* 2026. DOI: 10.64898/2026.04.07.716960

---

## Contact

**Meng Zhang**  
Department of Clinical Laboratory, Tongde Hospital of Zhejiang Province  
Email: zh3028@zuaa.zju.edu.cn  
ORCID: 0000-0002-2127-3126

---

## License

This project is licensed under the MIT License.
