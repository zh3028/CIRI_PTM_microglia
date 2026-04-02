# 01_preprocessing.R
# 数据预处理：读取、质控、SCTransform、降维、聚类
# 适用于 D1, D2, D3, D4 (D5 为 bulk，单独处理)
setwd("D:/GEO_data")
library(Seurat)
library(dplyr)
library(data.table)  # 用于 D5

# 创建目录
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ========== D1: GSE174574 ==========
untar("GSE174574_RAW.tar", exdir = "GSE174574_RAW")
# 假设已手动整理子文件夹（每个样本一个文件夹，内含三个标准文件）
sample_dirs <- list.dirs("GSE174574_RAW", recursive = FALSE, full.names = TRUE)
sample_ids <- basename(sample_dirs)
seurat_list <- list()
for(i in seq_along(sample_dirs)) {
  counts <- Read10X(sample_dirs[i])
  obj <- CreateSeuratObject(counts, project = sample_ids[i], min.cells = 3, min.features = 200)
  seurat_list[[sample_ids[i]]] <- obj
}
obj_D1 <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_ids, project = "GSE174574")
obj_D1$condition <- ifelse(obj_D1$orig.ident %in% c("GSM5319987","GSM5319988","GSM5319989"), "Sham", "MCAO")
obj_D1$condition <- factor(obj_D1$condition, levels = c("Sham", "MCAO"))
obj_D1[["percent.mt"]] <- PercentageFeatureSet(obj_D1, pattern = "^mt-")
obj_D1[["percent.rb"]] <- PercentageFeatureSet(obj_D1, pattern = "^Hb[ab]-")
obj_D1 <- subset(obj_D1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                   nCount_RNA > 500 & nCount_RNA < 30000 &
                   percent.mt < 20 & percent.rb < 5)
obj_D1 <- SCTransform(obj_D1, vars.to.regress = "percent.mt", verbose = FALSE)
obj_D1 <- RunPCA(obj_D1, npcs = 50)
obj_D1 <- RunUMAP(obj_D1, dims = 1:30, seed.use = 2026)
obj_D1 <- FindNeighbors(obj_D1, dims = 1:30)
obj_D1 <- FindClusters(obj_D1, resolution = 0.8, random.seed = 2026)
saveRDS(obj_D1, "data/processed/D1_raw_seurat.rds")

# ========== D2: GSE227651 ==========
untar("GSE227651_RAW.tar", exdir = "GSE227651_RAW")
sample_dirs <- list.dirs("GSE227651_RAW", recursive = FALSE, full.names = TRUE)
sample_ids <- basename(sample_dirs)
seurat_list <- list()
for(i in seq_along(sample_dirs)) {
  counts <- Read10X(sample_dirs[i])
  obj <- CreateSeuratObject(counts, project = sample_ids[i], min.cells = 3, min.features = 200)
  seurat_list[[sample_ids[i]]] <- obj
}
obj_D2 <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_ids, project = "GSE227651")
timepoint_map <- c("GSM7104632" = "d1", "GSM7104633" = "d3", "GSM7104634" = "d7", "GSM7104635" = "Sham")
obj_D2$timepoint <- factor(timepoint_map[obj_D2$orig.ident], levels = c("Sham", "d1", "d3", "d7"))
obj_D2[["percent.mt"]] <- PercentageFeatureSet(obj_D2, pattern = "^mt-")
obj_D2[["percent.rb"]] <- PercentageFeatureSet(obj_D2, pattern = "^Hb[ab]-")
obj_D2 <- subset(obj_D2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                   nCount_RNA > 500 & nCount_RNA < 30000 &
                   percent.mt < 20 & percent.rb < 5)
obj_D2 <- SCTransform(obj_D2, vars.to.regress = "percent.mt", verbose = FALSE)
obj_D2 <- RunPCA(obj_D2, npcs = 50)
obj_D2 <- RunUMAP(obj_D2, dims = 1:30, seed.use = 2026)
obj_D2 <- FindNeighbors(obj_D2, dims = 1:30)
obj_D2 <- FindClusters(obj_D2, resolution = 0.8, random.seed = 2026)
saveRDS(obj_D2, "data/processed/D2_raw_seurat.rds")

# ========== D3: GSE245386 ==========
untar("GSE245386_RAW.tar", exdir = "GSE245386_RAW")
# 读取所有样本，不预筛选（因为文件夹名不含 WTC/WTM，后续根据 orig.ident 分组）
sample_dirs <- list.dirs("GSE245386_RAW", recursive = FALSE, full.names = TRUE)
seurat_list <- list()
for(i in seq_along(sample_dirs)) {
  counts <- Read10X(sample_dirs[i])
  obj <- CreateSeuratObject(counts, project = basename(sample_dirs[i]), min.cells = 3, min.features = 200)
  seurat_list[[basename(sample_dirs[i])]] <- obj
}
obj_D3 <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "GSE245386")
# 根据 orig.ident 中的 "WTC" 或 "WTM" 设置 condition
obj_D3$condition <- ifelse(grepl("WTC", obj_D3$orig.ident), "Sham", "MCAO")
obj_D3$condition <- factor(obj_D3$condition, levels = c("Sham", "MCAO"))
obj_D3[["percent.mt"]] <- PercentageFeatureSet(obj_D3, pattern = "^mt-")
obj_D3[["percent.rb"]] <- PercentageFeatureSet(obj_D3, pattern = "^Hb[ab]-")
obj_D3 <- subset(obj_D3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                   nCount_RNA > 500 & nCount_RNA < 30000 &
                   percent.mt < 20 & percent.rb < 5)
obj_D3 <- SCTransform(obj_D3, vars.to.regress = "percent.mt", verbose = FALSE)
obj_D3 <- RunPCA(obj_D3, npcs = 50)
obj_D3 <- RunUMAP(obj_D3, dims = 1:30, seed.use = 2026)
obj_D3 <- FindNeighbors(obj_D3, dims = 1:30)
obj_D3 <- FindClusters(obj_D3, resolution = 0.8, random.seed = 2026)
saveRDS(obj_D3, "data/processed/D3_raw_seurat.rds")

# ========== D4: GSE267240 ==========
untar("GSE267240_RAW.tar", exdir = "GSE267240_RAW")
sample_dirs <- list.dirs("GSE267240_RAW", recursive = FALSE, full.names = TRUE)
sample_ids <- basename(sample_dirs)
seurat_list <- list()
for(i in seq_along(sample_dirs)) {
  counts <- Read10X(sample_dirs[i])
  obj <- CreateSeuratObject(counts, project = sample_ids[i], min.cells = 3, min.features = 200)
  seurat_list[[sample_ids[i]]] <- obj
}
obj_D4 <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_ids, project = "GSE267240")
obj_D4$sex <- ifelse(grepl("F_", obj_D4$orig.ident), "F", "M")
obj_D4$condition <- ifelse(grepl("Stroke", obj_D4$orig.ident), "Stroke", "Sham")
obj_D4$condition <- factor(obj_D4$condition, levels = c("Sham", "Stroke"))
obj_D4[["percent.mt"]] <- PercentageFeatureSet(obj_D4, pattern = "^mt-")
obj_D4 <- subset(obj_D4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                   nCount_RNA > 500 & nCount_RNA < 30000 &
                   percent.mt < 20)
obj_D4 <- SCTransform(obj_D4, vars.to.regress = "percent.mt", verbose = FALSE)
obj_D4 <- RunPCA(obj_D4, npcs = 50)
obj_D4 <- RunUMAP(obj_D4, dims = 1:30, seed.use = 2026)
obj_D4 <- FindNeighbors(obj_D4, dims = 1:30)
obj_D4 <- FindClusters(obj_D4, resolution = 0.8, random.seed = 2026)
saveRDS(obj_D4, "data/processed/D4_raw_seurat.rds")

# ========== D5: GSE319237 (bulk) ==========
if(file.exists("GSE319237_bulkRNA_timeseries_Day6Igf1_gene_counts.txt.gz")) {
  d5 <- fread("GSE319237_bulkRNA_timeseries_Day6Igf1_gene_counts.txt.gz", data.table = FALSE)
  rownames(d5) <- d5[,1]
  d5 <- d5[,-1]
  saveRDS(d5, "data/processed/D5_bulk_expression.rds")
} else {
  warning("D5 file not found, please download manually.")
}