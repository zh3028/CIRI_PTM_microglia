# 05_validation.R
# D3 独立验证
setwd("D:/GEO_data")
library(Seurat)

W_D1 <- as.matrix(read.csv("results/nmf_results/D1_W_matrix_correct.csv", row.names = 1))
ptm_genes <- readLines("results/tables/ptm_genes_final.txt")
obj_D3 <- readRDS("data/processed/D3_microglia.rds")

expr <- GetAssayData(obj_D3, assay = "SCT", layer = "data")[ptm_genes, , drop = FALSE]
X <- t(expr)
X <- X[rowSums(X) > 0, ]
common <- intersect(rownames(W_D1), colnames(X))
W_aligned <- W_D1[common, , drop = FALSE]
X_aligned <- X[, common, drop = FALSE]
H <- X_aligned %*% W_aligned %*% solve(t(W_aligned) %*% W_aligned)
H[H < 0] <- 0
module_assign <- apply(H, 1, which.max)
obj_D3$module_projected <- factor(module_assign, levels = 1:3, labels = paste0("M",1:3))
saveRDS(obj_D3, "data/processed/D3_microglia_with_projection.rds")

# D1 MCAO 组模块比例
obj_D1 <- readRDS("data/processed/D1_microglia_with_modules.rds")
d1_mcao <- obj_D1$module[obj_D1$condition == "MCAO" & !is.na(obj_D1$module)]
d1_prop <- table(d1_mcao) / length(d1_mcao)
d3_prop <- table(obj_D3$module_projected) / ncol(obj_D3)
cos_sim <- sum(d3_prop * d1_prop) / sqrt(sum(d3_prop^2) * sum(d1_prop^2))

# 保存结果
val_res <- data.frame(
  dataset = "D3",
  cosine_similarity = cos_sim,
  prop_M1_D3 = d3_prop["M1"], prop_M2_D3 = d3_prop["M2"], prop_M3_D3 = d3_prop["M3"],
  prop_M1_D1 = d1_prop["1"], prop_M2_D1 = d1_prop["2"], prop_M3_D1 = d1_prop["3"]
)
write.csv(val_res, "results/tables/D3_validation_results.csv", row.names = FALSE)

# 条形图
prop_data <- data.frame(
  dataset = c(rep("D1 MCAO",3), rep("D3",3)),
  module = rep(c("M1","M2","M3"),2),
  proportion = c(d1_prop["1"], d1_prop["2"], d1_prop["3"], d3_prop["M1"], d3_prop["M2"], d3_prop["M3"])
)
library(ggplot2)
pdf("results/figures/D3_validation_barplot.pdf", width = 5, height = 4)
ggplot(prop_data, aes(x = module, y = proportion, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("D1 MCAO" = "#4DAF4A", "D3" = "#E41A1C")) +
  labs(x = "Module", y = "Proportion of cells") + theme_minimal() +
  theme(legend.title = element_blank())
dev.off()