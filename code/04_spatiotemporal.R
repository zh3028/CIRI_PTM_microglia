# 04_spatiotemporal.R
# D2 投影、模块比例趋势、得分动态、卡方检验
setwd("D:/GEO_data")
library(Seurat)
library(dplyr)
library(ggplot2)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# 加载 D1 权重矩阵和 D2 小胶质细胞
W_D1 <- as.matrix(read.csv("results/nmf_results/D1_W_matrix_correct.csv", row.names = 1))
obj_D2 <- readRDS("data/processed/D2_microglia_loose.rds")
ptm_genes <- readLines("results/tables/ptm_genes_final.txt")

# 提取表达矩阵并投影
expr <- GetAssayData(obj_D2, assay = "SCT", layer = "data")[ptm_genes, , drop = FALSE]
X <- t(expr)
X <- X[rowSums(X) > 0, ]
common <- intersect(rownames(W_D1), colnames(X))
W_aligned <- W_D1[common, , drop = FALSE]
X_aligned <- X[, common, drop = FALSE]
H <- X_aligned %*% W_aligned %*% solve(t(W_aligned) %*% W_aligned)
H[H < 0] <- 0
colnames(H) <- paste0("M",1:3)
rownames(H) <- rownames(X_aligned)
obj_D2 <- obj_D2[, rownames(H)]
obj_D2$module_projected <- factor(apply(H, 1, which.max), levels = 1:3, labels = paste0("M",1:3))
saveRDS(obj_D2, "data/processed/D2_microglia_loose_with_projection.rds")

# 模块比例
prop <- obj_D2@meta.data %>% group_by(timepoint, module_projected) %>%
  summarise(count = n()) %>% group_by(timepoint) %>% mutate(proportion = count / sum(count))
write.csv(prop, "results/tables/D2_module_proportions.csv", row.names = FALSE)

# 绘图：比例趋势
pdf("results/figures/D2_module_proportion_trend.pdf", width = 8, height = 5)
ggplot(prop, aes(x = timepoint, y = proportion, color = module_projected, group = module_projected)) +
  geom_line(size = 1) + geom_point(size = 2) +
  labs(x = "Time after tMCAO", y = "Cell proportion", color = "Module") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 卡方检验
chisq_result <- chisq.test(table(obj_D2$module_projected, obj_D2$timepoint))
cat("Chi-square test: X-squared =", chisq_result$statistic, ", df =", chisq_result$parameter, ", p =", chisq_result$p.value, "\n")

# 模块得分动态（按时间点平均得分）
score_by_time <- data.frame(
  timepoint = obj_D2$timepoint,
  M1 = H[, "M1"],
  M2 = H[, "M2"],
  M3 = H[, "M3"]
) %>% group_by(timepoint) %>% summarise(mean_M1 = mean(M1), mean_M2 = mean(M2), mean_M3 = mean(M3))
write.csv(score_by_time, "results/tables/D2_module_scores_by_timepoint.csv", row.names = FALSE)

# 绘图：得分趋势
score_long <- data.frame(
  timepoint = rep(score_by_time$timepoint, 3),
  module = rep(c("M1","M2","M3"), each = nrow(score_by_time)),
  score = c(score_by_time$mean_M1, score_by_time$mean_M2, score_by_time$mean_M3)
)
pdf("results/figures/D2_module_scores_by_timepoint.pdf", width = 6, height = 4)
ggplot(score_long, aes(x = timepoint, y = score, color = module, group = module)) +
  geom_line(size = 1) + geom_point(size = 2) +
  labs(x = "Time after tMCAO", y = "Mean module score", color = "Module") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()