# PTM enzyme co-expression modules in microglia after cerebral ischemia-reperfusion injury

[![R](https://img.shields.io/badge/R-4.3.3-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the computational pipeline for the manuscript:
> **A PTM Regulatory Enzyme Co-expression Code Defines Microglial Functional Heterogeneity in Cerebral Ischemia-Reperfusion Injury**

We identified three PTM enzyme co-expression modules (metabolic stress, pro-inflammatory, and reparative) that exhibit distinct spatiotemporal dynamics in microglia after transient middle cerebral artery occlusion (tMCAO). The analysis integrates five public GEO datasets.

---

## 📁 Repository structure
```
├── code/             # Analysis scripts (R)
│   ├── 01_preprocessing.R
│   ├── 02_microglia_filter.R
│   ├── 03_nmf_analysis.R
│   ├── 04_spatiotemporal.R
│   ├── 05_validation.R
│   ├── 06_sex_analysis.R
│   ├── 07_purity_validation.R
│   └── 08_figures.R
└── LICENSE
```
