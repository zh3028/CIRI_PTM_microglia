# 06_sex_analysis.R
# D4 性别差异分析
setwd("D:/GEO_data")
library(Seurat)
library(dplyr)
library(ggplot2)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

W_D1 <- as.matrix(read.csv("results/nmf_results/D1_W_matrix_correct.csv", row.names = 1))
ptm_genes <- readLines("results/tables/ptm_genes_final.txt")
obj_D4 <- readRDS("data/processed/D4_microglia.rds")

# 投影
expr <- GetAssayData(obj_D4, assay = "SCT", layer = "data")[ptm_genes, , drop = FALSE]
X <- t(expr)
X <- X[rowSums(X) > 0, ]
common <- intersect(rownames(W_D1), colnames(X))
W_aligned <- W_D1[common, , drop = FALSE]
X_aligned <- X[, common, drop = FALSE]
H <- X_aligned %*% W_aligned %*% solve(t(W_aligned) %*% W_aligned)
H[H < 0] <- 0
module_assign <- apply(H, 1, which.max)
obj_D4$module <- factor(module_assign, levels = 1:3, labels = paste0("M",1:3))
saveRDS(obj_D4, "data/processed/D4_with_projection.rds")

# 性别差异（仅 stroke 组）
stroke <- subset(obj_D4, condition == "Stroke")
tbl <- table(stroke$sex, stroke$module)
chisq <- chisq.test(tbl)
fisher <- fisher.test(tbl, simulate.p.value = TRUE, B = 10000)

# bootstrap 验证
set.seed(2026)
boot_p <- numeric(1000)
for(b in 1:1000) {
  fem <- sample(stroke$module[stroke$sex == "F"], replace = TRUE)
  mal <- sample(stroke$module[stroke$sex == "M"], replace = TRUE)
  bt <- table(c(fem, mal), c(rep("F",length(fem)), rep("M",length(mal))))
  if(nrow(bt)==3 && ncol(bt)==2) boot_p[b] <- chisq.test(bt)$p.value
}
boot_p <- na.omit(boot_p)
boot_median <- median(boot_p)
boot_ci <- quantile(boot_p, c(0.025,0.975))

# 保存结果
sex_res <- data.frame(
  test = c("Chi-square","Fisher exact","Bootstrap median","Bootstrap 2.5%","Bootstrap 97.5%"),
  p_value = c(chisq$p.value, fisher$p.value, boot_median, boot_ci[1], boot_ci[2])
)
write.csv(sex_res, "results/tables/D4_sex_analysis.csv", row.names = FALSE)

# 条形图
prop_df <- as.data.frame(tbl) %>% group_by(Var1) %>% mutate(proportion = Freq / sum(Freq))
pdf("results/figures/D4_sex_module_barplot.pdf", width = 6, height = 4)
ggplot(prop_df, aes(x = Var1, y = proportion, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sex (Stroke group)", y = "Proportion", fill = "Module") + theme_minimal()
dev.off()