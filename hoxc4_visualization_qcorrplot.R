# Hoxc4通信矩阵可视化 - 使用qcorrplot和linkET

cat("\n=== Hoxc4通信矩阵可视化测试 ===\n")

# 加载必要的依赖包
cat("\n加载必要的依赖包...\n")
required_packages <- c("ggplot2", "reshape2", "vegan", "RColorBrewer", "linkET", "dplyr")

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

# 创建Mantel测试结果数据框
mantel <- data.frame(
  x = rep("DOCK5", n_cells),
  y = cell_types,
  pd = cut(mantel_result$signif, breaks = c(0, 0.01, 0.05, 1), labels = c("< 0.01", "0.01-0.05", "> 0.05")),
  rd = cut(mantel_result$statistic, breaks = c(0, 0.2, 0.4, 1), labels = c("< 0.2", "0.2-0.4", "> 0.4"))
)

# 确保mantel数据框有足够的行
if (nrow(mantel) < n_cells) {
  # 为缺失的细胞类型生成随机的Mantel测试结果
  for (i in (nrow(mantel) + 1):n_cells) {
    random_p <- runif(1, 0.001, 0.1)
    random_r <- runif(1, 0.3, 0.9)
    mantel <- rbind(mantel, data.frame(
      x = "DOCK5",
      y = cell_types[i],
      pd = cut(random_p, breaks = c(0, 0.01, 0.05, 1), labels = c("< 0.01", "0.01-0.05", "> 0.05")),
      rd = cut(random_r, breaks = c(0, 0.2, 0.4, 1), labels = c("< 0.2", "0.2-0.4", "> 0.4"))
    ))
  }
}

# 使用用户指定的代码绘制
cat("\n使用qcorrplot绘制可视化...\n")

# 定义nice_curvature函数（linkET包中可能需要）
nice_curvature <- function() {
  return(0.1)
}

# 定义color_pal函数
color_pal <- function(n) {
  return(RColorBrewer::brewer.pal(min(n, 8), "Set1"))
}

# 执行用户指定的绘图代码
p <- qcorrplot(correlate(sym_matrix_high), type = "upper", diag = FALSE) + 
  # 使用 geom_tile() 替代 geom_square()，或者提供完整的美学映射 
  geom_tile() +  # 替代方案1：使用 geom_tile() 
  # 或者提供完整的美学映射给 geom_square() 
  # geom_square(aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5)) +  # 替代方案2 
  
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) + 
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) + 
  scale_size_manual(values = c(0.5, 1, 2)) + 
  scale_colour_manual(values = color_pal(3)) + 
  guides(size = guide_legend(title = "Mantel's r", 
                             override.aes = list(colour = "grey35"), 
                             order = 2), 
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1), 
         fill = guide_colorbar(title = "Pearson's r", order = 3))

# 添加标题
p <- p + ggtitle("Hoxc4基因通信网络相关性分析")

# 保存图形
output_file <- file.path(output_dir, "hoxc4_qcorrplot_analysis.png")
ggsave(output_file, plot = p, width = 14, height = 12, dpi = 300)
cat("\n✓ 可视化已保存至:", output_file)

# 显示结果摘要
cat("\n\n===== 可视化完成 =====\n")
cat("生成的文件:\n")
cat("1.", output_file, "\n")
