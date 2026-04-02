# 02_microglia_filter.R
# 小胶质细胞双重纯度过滤
setwd("D:/GEO_data")
library(Seurat)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# 定义过滤函数
filter_microglia <- function(obj, pos_genes = c("P2ry12", "Tmem119", "Cx3cr1"),
                             neg_genes = c("Cd68", "Adgre1", "Ly6c")) {
  expr <- GetAssayData(obj, assay = "SCT", layer = "data")
  pos_expr <- expr[pos_genes, , drop = FALSE] > 0
  pos_cells <- colSums(pos_expr) >= 2
  neg_exist <- intersect(neg_genes, rownames(expr))
  if(length(neg_exist) > 0) {
    neg_expr <- expr[neg_exist, , drop = FALSE] > 0
    neg_cells <- colSums(neg_expr) == 0
  } else {
    neg_cells <- rep(TRUE, ncol(expr))
  }
  keep <- pos_cells & neg_cells
  obj_filtered <- subset(obj, cells = colnames(obj)[keep])
  return(obj_filtered)
}

# 对 D1 过滤
obj_D1 <- readRDS("data/processed/D1_raw_seurat.rds")
obj_D1_micro <- filter_microglia(obj_D1)
saveRDS(obj_D1_micro, "data/processed/D1_microglia.rds")

# 对 D2 过滤（宽松条件：至少1个阳性基因）
obj_D2 <- readRDS("data/processed/D2_raw_seurat.rds")
expr <- GetAssayData(obj_D2, assay = "SCT", layer = "data")
pos_genes <- c("P2ry12", "Tmem119", "Cx3cr1")
pos_cells <- colSums(expr[pos_genes, , drop = FALSE] > 0) >= 1
neg_genes_all <- c("Cd68", "Adgre1", "Ly6c")
neg_genes_exist <- intersect(neg_genes_all, rownames(expr))
if(length(neg_genes_exist) > 0) {
  neg_cells <- colSums(expr[neg_genes_exist, , drop = FALSE] > 0) == 0
} else {
  neg_cells <- rep(TRUE, ncol(expr))
}
keep <- pos_cells & neg_cells
obj_D2_micro <- subset(obj_D2, cells = colnames(obj_D2)[keep])
saveRDS(obj_D2_micro, "data/processed/D2_microglia_loose.rds")

# 对 D3 过滤
obj_D3 <- readRDS("data/processed/D3_raw_seurat.rds")
obj_D3_micro <- filter_microglia(obj_D3)
saveRDS(obj_D3_micro, "data/processed/D3_microglia.rds")

# D4 已分选，无需过滤，但为了统一命名，直接复制
obj_D4 <- readRDS("data/processed/D4_raw_seurat.rds")
saveRDS(obj_D4, "data/processed/D4_microglia.rds")