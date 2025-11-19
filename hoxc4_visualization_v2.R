# Hoxc4通信矩阵可视化简化脚本

cat("\n=== Hoxc4通信矩阵可视化测试 ===\n")

# 加载必要的依赖包
cat("\n加载必要的依赖包...\n")
required_packages <- c("ggplot2", "reshape2", "vegan", "RColorBrewer", "dplyr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat("正在安装", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 设置输出路径
output_dir <- file.path(getwd(), "hoxc4_visualization")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
cat("输出目录:", output_dir, "\n")

# 创建模拟数据
cat("\n创建模拟数据...\n")
set.seed(123)
n_cells <- 8
cell_types <- c("B细胞", "树突状细胞", "内皮细胞", "成纤维细胞", "肥大细胞", "髓系细胞", "卵巢癌细胞", "浆细胞")

# 创建模拟的通信矩阵
create_matrix <- function() {
  # 创建基础矩阵
  mat <- matrix(runif(n_cells*n_cells, 0, 0.3), nrow = n_cells, ncol = n_cells)
  # 添加一些强相关性
  for (i in 1:5) {
    for (j in 1:5) {
      mat[i, j] <- mat[i, j] + runif(1, 0.4, 0.7)
    }
  }
  # 使矩阵对称
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 0
  rownames(mat) <- colnames(mat) <- cell_types
  return(mat)
}

# 创建高低表达组矩阵
sym_matrix_high <- create_matrix()
sym_matrix_low <- 0.6 * sym_matrix_high

# 执行Mantel测试
cat("\n执行Mantel测试...\n")
mantel_result <- vegan::mantel(
  as.dist(1 - sym_matrix_high),
  as.dist(1 - sym_matrix_low),
  method = "pearson",
  permutations = 999
)

# 显示Mantel测试结果
cat("Mantel测试结果:")
cat("\n相关系数(r):", mantel_result$statistic)
cat("\np值:", mantel_result$signif)

# 创建相关矩阵
correlation_matrix <- cor(sym_matrix_high, method = "pearson")

# 创建热图
cat("\n创建热图可视化...\n")
corr_df <- reshape2::melt(correlation_matrix)
colnames(corr_df) <- c("Var1", "Var2", "correlation")

# 创建热图 - 参照用户提供的理想效果
p <- ggplot(corr_df, aes(x = Var2, y = Var1, fill = correlation))
p <- p + geom_tile(color = "white", size = 0.5)
p <- p + scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "Pearson相关系数",
                        limits = c(-1, 1))
p <- p + theme_minimal()
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
p <- p + ggtitle("Hoxc4基因通信网络相关性分析")
p <- p + labs(x = "", y = "")

# 添加相关性连线
sig_corr <- which(lower.tri(correlation_matrix), arr.ind = TRUE)
if (nrow(sig_corr) > 0) {
  link_data <- data.frame(
    x = colnames(correlation_matrix)[sig_corr[, 2]],
    y = rownames(correlation_matrix)[sig_corr[, 1]],
    corr = correlation_matrix[sig_corr]
  )
  
  p <- p + geom_segment(data = link_data,
                aes(x = x, y = y, xend = y, yend = x, 
                    color = corr, linewidth = abs(corr)),
                alpha = 0.8,
                inherit.aes = FALSE)
  p <- p + scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, guide = guide_colorbar(title = "相关系数"))
  p <- p + scale_linewidth_continuous(range = c(0.3, 1.5), guide = guide_legend(title = "强度"))
}

# 保存热图
output_file <- file.path(output_dir, "hoxc4_correlation_analysis.png")
ggsave(output_file, plot = p, width = 12, height = 10, dpi = 300)
cat("\n✓ 热图已保存至:", output_file)

# 创建增强版对比图 - 参照用户期望的效果
cat("\n\n创建增强版对比图...\n")

# 计算每个细胞类型的总体通信强度
high_mean_intensity <- apply(sym_matrix_high, 1, mean)
low_mean_intensity <- apply(sym_matrix_low, 1, mean)

# 计算差异
diff_intensity <- high_mean_intensity - low_mean_intensity

# 创建差异数据框
diff_df <- data.frame(
  cell_type = cell_types,
  intensity = diff_intensity
)

# 添加模拟的Mantel测试结果（为了显示效果）
mantel_p_values <- runif(n_cells, 0.001, 0.1)
mantel_p_groups <- ifelse(mantel_p_values < 0.01, "< 0.01", 
                         ifelse(mantel_p_values < 0.05, "0.01-0.05", "> 0.05"))

mantel_r_values <- runif(n_cells, 0.3, 0.9)
mantel_r_groups <- ifelse(mantel_r_values <= 0.4, "0.2-0.4", "> 0.4")

diff_df$p_value <- mantel_p_values
diff_df$p_group <- mantel_p_groups
diff_df$r_value <- mantel_r_values
diff_df$r_group <- mantel_r_groups

# 创建增强版对比图
# 主图：柱状图显示通信强度差异
p_compare <- ggplot(diff_df, aes(x = cell_type, y = intensity, fill = intensity))
p_compare <- p_compare + geom_bar(stat = "identity")

# 添加颜色映射
p_compare <- p_compare + scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Pearson相关系数")

# 添加Mantel测试结果点
p_compare <- p_compare + geom_point(data = diff_df, 
                        aes(x = cell_type, y = 0, color = p_group),
                        size = 5, position = position_dodge(width = 0.9))

# 添加模拟的连线效果
# 模拟DOCK5与细胞类型的连线
dock5_connections <- data.frame(
  x = rep("DOCK5", n_cells),
  y = rep(1, n_cells),
  xend = cell_types,
  yend = diff_intensity / 2,
  strength = mantel_r_values,
  p_group = mantel_p_groups
)

# 添加连线
p_compare <- p_compare + geom_segment(data = dock5_connections,
                         aes(x = x, y = y, xend = xend, yend = yend, 
                             color = p_group, linewidth = strength),
                         alpha = 0.7, inherit.aes = FALSE)

# 添加图例
p_compare <- p_compare + scale_color_manual(name = "Mantel检验p值",
                         values = c("< 0.01" = "red", "0.01-0.05" = "orange", "> 0.05" = "gray"),
                         breaks = c("< 0.01", "0.01-0.05", "> 0.05"))

p_compare <- p_compare + scale_linewidth_continuous(name = "Mantel检验r值",
                             breaks = c(0.4, 0.6, 0.8),
                             labels = c("0.2-0.4", "0.4-0.6", "> 0.8"))

# 设置主题
p_compare <- p_compare + theme_minimal()
p_compare <- p_compare + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# 添加标题和标签
p_compare <- p_compare + ggtitle("Hoxc4基因高/低表达组通信网络对比")
p_compare <- p_compare + labs(x = "细胞类型", y = "通信强度差异")

# 保存对比图
output_compare_file <- file.path(output_dir, "hoxc4_expression_comparison.png")
ggsave(output_compare_file, plot = p_compare, width = 14, height = 6, dpi = 300)
cat("\n✓ 对比图已保存至:", output_compare_file)

# 显示结果摘要
cat("\n\n===== 可视化完成 =====\n")
cat("生成的文件:\n")
cat("1.", output_file, "\n")
cat("2.", output_compare_file, "\n")
