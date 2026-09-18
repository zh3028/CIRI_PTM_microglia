# Data and Code for Two Related Studies on Diabetic Retinopathy

## Overview

This repository contains analysis scripts for two related studies on diabetic retinopathy (DR):

- **Study 1 (in submission):** Integrative single-cell and machine learning analysis identifies a ferroptosis-related three-gene diagnostic signature for diabetic retinopathy
- **Study 2 (in preparation):** Therapeutic mechanism of *Tetrastigma hemsleyanum* in diabetic retinopathy: targeting HMOX1/CTSD-mediated ferroptosis in pathogenic macrophages

## ScienceDB DOI

All processed data (including Seurat objects, machine learning models, expression matrices, and ROC objects) are available at ScienceDB:
https://doi.org/10.57760/sciencedb.0134m

## Repository Contents

- `01_scripts/`: R scripts used for all analyses
- `05_session_info/`: R session information for reproducibility
- `README.md`: This file

## Data Sources

- GSE165784: Human PDR fibrovascular membranes (scRNA-seq)
- GSE209872: Rat STZ-induced DR retina (scRNA-seq)
- GSE102485: Human PDR bulk RNA-seq (training cohort)
- GSE160306: Human DR bulk RNA-seq (external validation)
- GSE60436: Human PDR microarray (external validation)

## Key Findings

### Study 1

1. Identified 11 macrophage subpopulations in human PDR
2. SPP1+/HMOX1+ macrophages exhibit ferroptosis activation (HMOX1, CTSD, FTH1)
3. Three-gene diagnostic signature (FTH1, LDHA, CTSZ) with 10-fold CV AUC=0.800 (training), 0.889 (GSE60436 after Z-score), 0.705 (FTH1 alone in GSE160306)
4. FTH1 correlates with monocyte/macrophage infiltration (r=0.316, P=0.015)

### Study 2

1. CellChat analysis reveals SPP1_Mac communicates with pericytes (0.907), fibroblasts (0.819), and endothelial cells (0.658)
2. Molecular docking shows active compounds of *Tetrastigma hemsleyanum* bind to HMOX1, CTSD, and SPP1
3. Cross-species conservation of pathogenic signatures in rat STZ-DR model
4. Wet-lab validation (in progress): HPLC quantification of ocular penetration, cell and animal experiments

## Contact

Meng Zhang, Department of Clinical Laboratory, Tongde Hospital of Zhejiang Province
Email: zh3028@zuaa.zju.edu.cn
