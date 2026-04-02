# 03_nmf_analysis.R
# NMF 模块识别、稳定性验证、功能注释
setwd("D:/GEO_data")
library(Seurat)
library(NMF)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(clue)          # 用于 solve_LSAP

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/nmf_results", recursive = TRUE, showWarnings = FALSE)

# 加载 D1 小胶质细胞
obj_D1 <- readRDS("data/processed/D1_microglia.rds")

# PTM 基因列表（25个，过滤后保留22个）
ptm_genes <- c("Acod1","Gsr","Gclc","Gclm","Suclg1","Sucla2","Sdha","Sdhb","Kat2a","Kat2b","Sirt5",
               "Ldha","Ldhb","Slc16a1","Slc16a3","Ep300","Hdac1","Hdac2","Hdac3","Sirt1","Sirt2","Sirt3","Padi4")
expr_data <- GetAssayData(obj_D1, assay = "SCT", layer = "data")
available <- intersect(ptm_genes, rownames(expr_data))
pct_exp <- rowSums(expr_data[available, ] > 0) / ncol(expr_data) * 100
ptm_final <- available[pct_exp >= 5]
writeLines(ptm_final, "results/tables/ptm_genes_final.txt")

# 准备 NMF 输入矩阵（细胞 × 基因）
X <- t(GetAssayData(obj_D1, assay = "SCT", layer = "data")[ptm_final, ])
X <- X[rowSums(X) > 0, ]  # 移除全零细胞
X_dense <- as.matrix(X)
write.csv(X_dense, "results/nmf_results/D1_nmf_input_matrix.csv", row.names = TRUE)

# nrun 稳定性测试
set.seed(2026)
nrun_test <- c(20,50,100)
W_list <- list()
for(n in nrun_test) {
  res <- nmf(X_dense, rank = 3, method = "brunet", nrun = n, seed = 2026)
  W_list[[as.character(n)]] <- basis(res)
}
cor_20_50 <- cor(as.vector(W_list[["20"]]), as.vector(W_list[["50"]]))
cor_20_100 <- cor(as.vector(W_list[["20"]]), as.vector(W_list[["100"]]))
cor_50_100 <- cor(as.vector(W_list[["50"]]), as.vector(W_list[["100"]]))
nrun_stab <- data.frame(comparison = c("20 vs 50","20 vs 100","50 vs 100"),
                        correlation = c(cor_20_50, cor_20_100, cor_50_100))
write.csv(nrun_stab, "results/nmf_results/nrun_stability_report.csv", row.names = FALSE)

# 秩选择（k=2:10）
res_rank <- nmf(X_dense, rank = 2:10, method = "brunet", nrun = 50, seed = 2026)
pdf("results/figures/D1_nmf_rank_selection.pdf", width = 10, height = 8)
plot(res_rank)
dev.off()
rank_summary <- summary(res_rank)
rank_eval <- data.frame(rank = 2:10, cophenetic = rank_summary$cophenetic, rss = rank_summary$rss)
write.csv(rank_eval, "results/nmf_results/rank_selection_metrics.csv", row.names = FALSE)

# 最终 NMF (k=3)
res_final <- nmf(X_dense, rank = 3, method = "brunet", nrun = 100, seed = 2026)
saveRDS(res_final, "results/nmf_results/D1_nmf_full.rds")
W_D1 <- basis(res_final)   # 基因 × 模块
H_D1 <- coef(res_final)    # 模块 × 细胞
write.csv(W_D1, "results/nmf_results/D1_W_matrix_correct.csv", row.names = TRUE)
write.csv(H_D1, "results/nmf_results/D1_H_matrix_correct.csv", row.names = TRUE)

# 细胞归属
module_assign <- apply(H_D1, 2, which.max)
names(module_assign) <- colnames(X_dense)
obj_D1$module <- NA
common <- intersect(names(module_assign), colnames(obj_D1))
obj_D1$module[common] <- factor(module_assign[common], levels = 1:3, labels = paste0("M",1:3))
saveRDS(obj_D1, "data/processed/D1_microglia_with_modules.rds")

# 交叉验证（5折）
set.seed(2026)
n_cells <- nrow(X_dense)
folds <- sample(rep(1:5, length.out = n_cells))
consistency <- numeric(5)
for(i in 1:5) {
  train_idx <- which(folds != i)
  test_idx <- which(folds == i)
  X_train <- X_dense[train_idx, ]
  X_test <- X_dense[test_idx, ]
  res_train <- nmf(X_train, rank = 3, method = "brunet", nrun = 50, seed = 2026)
  W_train <- basis(res_train)
  H_test <- t(solve(t(W_train) %*% W_train) %*% t(W_train) %*% t(X_test))
  H_test[H_test < 0] <- 0
  pred <- apply(H_test, 1, which.max)
  true <- module_assign[test_idx]
  confusion <- table(pred, true)
  matching <- solve_LSAP(-confusion, maximum = TRUE)
  pred_matched <- factor(pred, levels = 1:3, labels = as.character(1:3)[as.numeric(matching)])
  true_factor <- factor(true, levels = 1:3)
  consistency[i] <- sum(as.numeric(pred_matched) == as.numeric(true_factor)) / length(test_idx)
}
cv_report <- data.frame(fold = 1:5, consistency = consistency)
write.csv(cv_report, "results/nmf_results/cross_validation_report.csv", row.names = FALSE)

# 样本级验证（三个样本）
dir.create("results/nmf_results/D1_sample_level", recursive = TRUE, showWarnings = FALSE)
sample_names <- unique(obj_D1$orig.ident)
sample_cors <- numeric(length(sample_names))
for(s in seq_along(sample_names)) {
  cells <- WhichCells(obj_D1, expression = orig.ident == sample_names[s])
  mat_s <- X_dense[rownames(X_dense) %in% cells, , drop = FALSE]
  if(nrow(mat_s) < 50) next
  res_s <- nmf(mat_s, rank = 3, method = "brunet", nrun = 50, seed = 2026)
  W_s <- basis(res_s)
  cor_matrix_s <- matrix(0,3,3)
  for(i in 1:3) for(j in 1:3) cor_matrix_s[i,j] <- cor(W_s[,i], W_D1[,j])
  matching_s <- solve_LSAP(-cor_matrix_s, maximum = TRUE)
  W_s_matched <- W_s[, as.numeric(matching_s)]
  sample_cors[s] <- cor(as.vector(W_s_matched), as.vector(W_D1))
  write.csv(W_s, paste0("results/nmf_results/D1_sample_level/sample", s, "_W.csv"), row.names = TRUE)
}
sample_report <- data.frame(sample = sample_names, correlation = sample_cors)
write.csv(sample_report, "results/nmf_results/sample_level_validation.csv", row.names = FALSE)

# 共识图
pdf("results/figures/D1_nmf_consensusmap.pdf", width = 10, height = 8)
consensusmap(res_final, annCol = data.frame(module = obj_D1$module, row.names = colnames(X_dense)),
             main = "Consensus map of NMF modules")
dev.off()

# 差异表达与GO富集
obj_D1$module <- factor(obj_D1$module, levels = paste0("M",1:3))
markers <- FindAllMarkers(obj_D1, group.by = "module", logfc.threshold = 0.5, min.pct = 0.25,
                          only.pos = TRUE, random.seed = 2026)
write.csv(markers, "results/tables/D1_module_markers.csv", row.names = FALSE)
top_markers <- markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
write.csv(top_markers, "results/tables/D1_top_markers.csv", row.names = FALSE)

# GO 富集
convert_to_entrez <- function(gene_symbols) {
  map <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = TRUE)
  if(is.null(map)) return(character(0))
  return(map$ENTREZID)
}
for(mod in c("M1","M2","M3")) {
  genes <- markers %>% filter(cluster == mod, p_val_adj < 0.05) %>% pull(gene)
  if(length(genes) < 3) next
  entrez <- convert_to_entrez(genes)
  if(length(entrez) < 3) next
  ego <- enrichGO(gene = entrez, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  if(!is.null(ego) && nrow(ego) > 0) write.csv(as.data.frame(ego), paste0("results/tables/D1_module_", mod, "_GO.csv"), row.names = FALSE)
}