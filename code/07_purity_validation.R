# 07_purity_validation.R
# D5 纯度验证
setwd("D:/GEO_data")
library(Seurat)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

ptm_genes <- readLines("results/tables/ptm_genes_final.txt")

# 读取 D5 bulk 数据
if(file.exists("data/processed/D5_bulk_expression.rds")) {
  d5 <- readRDS("data/processed/D5_bulk_expression.rds")
  # 根据实际列名调整
  sham_col <- "bulkRNA_Sham_WT"
  stroke_col <- "bulkRNA_Day01_WT"
  common_genes <- intersect(ptm_genes, rownames(d5))
  d5_ptm <- d5[common_genes, c(sham_col, stroke_col)]
  d5_log2FC <- log2(d5_ptm[, stroke_col] + 1) - log2(d5_ptm[, sham_col] + 1)
} else {
  stop("D5 file not found")
}

# D1 中对应基因的 log2FC
obj_D1 <- readRDS("data/processed/D1_microglia.rds")
expr_D1 <- GetAssayData(obj_D1, assay = "SCT", layer = "data")[common_genes, ]
sham_cells <- WhichCells(obj_D1, expression = condition == "Sham")
stroke_cells <- WhichCells(obj_D1, expression = condition == "MCAO")
d1_log2FC <- log2(rowMeans(expr_D1[, stroke_cells], na.rm = TRUE) + 1) -
  log2(rowMeans(expr_D1[, sham_cells], na.rm = TRUE) + 1)
names(d1_log2FC) <- common_genes

# 相关性
cor_test <- cor.test(d1_log2FC, d5_log2FC, method = "spearman")
cor_report <- data.frame(method = "Spearman", correlation = cor_test$estimate, p_value = cor_test$p.value)
write.csv(cor_report, "results/tables/D1_D5_correlation.csv", row.names = FALSE)

# 散点图
pdf("results/figures/D1_D5_scatter.pdf", width = 6, height = 6)
plot(d1_log2FC, d5_log2FC, xlab = "D1 log2FC (scRNA-seq microglia)",
     ylab = "D5 log2FC (FACS-sorted microglia bulk)",
     main = paste("Spearman r =", round(cor_test$estimate, 3), "p =", format(cor_test$p.value, digits = 3)))
abline(0, 1, col = "red", lty = 2)
text(d1_log2FC, d5_log2FC, labels = names(d1_log2FC), pos = 3, cex = 0.5)
dev.off()