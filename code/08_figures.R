# 08_figures.R
# 生成所有主图和补充图（假设之前脚本已生成中间数据）
setwd("D:/GEO_data")
library(Seurat)
library(ggplot2)
library(pheatmap)
library(tidyr)        # 用于 pivot_longer
library(dplyr)

dir.create("results/figures/QC", recursive = TRUE, showWarnings = FALSE)

# 图 2A: D1 UMAP by condition
obj_D1 <- readRDS("data/processed/D1_microglia.rds")
pdf("results/figures/D1_umap_condition.pdf", width = 6, height = 5)
DimPlot(obj_D1, reduction = "umap", group.by = "condition", label = FALSE) +
  ggtitle("D1 microglia: Sham vs MCAO") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
dev.off()

# 图 2B: 标记基因热图
expr <- GetAssayData(obj_D1, assay = "SCT", layer = "data")
marker_genes <- c("P2ry12","Tmem119","Cx3cr1","Cd68","Adgre1","Ly6c")
marker_genes <- marker_genes[marker_genes %in% rownames(expr)]
pdf("results/figures/D1_marker_genes_heatmap.pdf", width = 8, height = 6)
DoHeatmap(obj_D1, features = marker_genes, group.by = "condition", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("Microglia purity validation")
dev.off()

# 图 3B: 模块权重热图
W <- as.matrix(read.csv("results/nmf_results/D1_W_matrix_correct.csv", row.names = 1))
rownames(W) <- gsub("^X", "", rownames(W))
colnames(W) <- c("M1","M2","M3")
W_norm <- t(scale(t(W)))
pdf("results/figures/D1_module_weight_heatmap.pdf", width = 6, height = 8)
pheatmap(W_norm, cluster_rows = TRUE, cluster_cols = FALSE,
         main = "PTM enzyme weights in modules",
         color = colorRampPalette(c("navy","white","firebrick3"))(50))
dev.off()

# 图 S1: 质控小提琴图（仅 D1-D3，D4 因缺少过滤对象跳过）
datasets <- list(
  D1 = list(raw = "data/processed/D1_raw_seurat.rds", filt = "data/processed/D1_microglia.rds"),
  D2 = list(raw = "data/processed/D2_raw_seurat.rds", filt = "data/processed/D2_microglia_loose.rds"),
  D3 = list(raw = "data/processed/D3_raw_seurat.rds", filt = "data/processed/D3_microglia.rds")
)
for(ds in names(datasets)) {
  raw <- readRDS(datasets[[ds]]$raw)
  filt <- readRDS(datasets[[ds]]$filt)
  raw_meta <- raw@meta.data; raw_meta$status <- "Before QC"
  filt_meta <- filt@meta.data; filt_meta$status <- "After QC"
  meta <- rbind(raw_meta[,c("nFeature_RNA","nCount_RNA","percent.mt","status")],
                filt_meta[,c("nFeature_RNA","nCount_RNA","percent.mt","status")])
  meta_long <- pivot_longer(meta, cols = c(nFeature_RNA, nCount_RNA, percent.mt),
                            names_to = "metric", values_to = "value")
  p <- ggplot(meta_long, aes(x = status, y = value, fill = status)) +
    geom_violin(scale = "width", trim = TRUE) + facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    labs(title = ds) + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(paste0("results/figures/QC/", ds, "_qc_violin.pdf"), p, width = 6, height = 8)
}

# 其他图已在相应脚本中生成，此处不再重复