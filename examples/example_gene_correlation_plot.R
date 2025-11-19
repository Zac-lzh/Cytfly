# example_gene_correlation_plot.R
# 示例脚本：演示如何使用plot_gene_specific_correlation函数

# 加载必要的包
library(Seurat)

# 加载CellChatCompare包中的自定义函数
source("../R/plot_gene_specific_correlation.R")

# 示例1：使用Seurat对象和指定基因生成相关性热图
tryCatch({
  # 假设你已经有一个名为seurat_obj的Seurat对象
  # 如果你需要创建一个测试对象，取消下面的注释
  
  # 创建一个小型测试Seurat对象
  # 创建随机表达矩阵
  n_cells <- 100
  n_genes <- 200
  set.seed(123)
  
  # 确保HOXC4基因在矩阵中
  gene_names <- c("HOXC4", paste0("Gene_", 1:(n_genes-1)))
  
  # 创建基础counts矩阵
  counts_matrix <- matrix(rpois(n_cells * n_genes, lambda = 3), nrow = n_genes)
  rownames(counts_matrix) <- gene_names
  colnames(counts_matrix) <- paste0("Cell_", 1:n_cells)
  
  # 为HOXC4基因设置差异表达
  hoxc4_counts <- rpois(n_cells, lambda = 5)
  # 使一部分细胞高表达HOXC4
  high_expr_indices <- sample(1:n_cells, 20)
  hoxc4_counts[high_expr_indices] <- hoxc4_counts[high_expr_indices] + 10
  counts_matrix["HOXC4", ] <- hoxc4_counts
  
  # 创建Seurat对象 - 使用公共API
  seurat_obj <- CreateSeuratObject(counts = counts_matrix)
  
  # 添加一些细胞类型信息 - 使用公共API
  seurat_obj <- AddMetaData(
    object = seurat_obj,
    metadata = sample(c("B细胞", "T细胞", "树突状细胞", "内皮细胞", "成纤维细胞"), 
                     size = n_cells, replace = TRUE),
    col.name = "cell_type"
  )
  
  # 运行基本的标准化和归一化
  seurat_obj <- NormalizeData(seurat_obj)
  
  # 确认HOXC4基因在对象中
  if (!"HOXC4" %in% rownames(seurat_obj)) {
    warning("HOXC4基因未在Seurat对象中找到，将使用模拟数据")
  } else {
    cat("HOXC4基因已成功添加到Seurat对象中\n")
  }
  
  # 运行基因特异性相关性分析
  cat("运行HOXC4基因的相关性分析...\n")
  results <- plot_gene_specific_correlation(
    seurat_object = seurat_obj,
    gene_name = "HOXC4",
    output_dir = "../output/hoxc4_correlation",
    top_percentile = 20,  # 前20%作为高表达组
    bottom_percentile = 20  # 后20%作为低表达组
  )
  
  cat("\n分析完成！结果存储在:", results$output_dir, "\n")
  cat("生成的文件包括:\n")
  cat("1. HOXC4_correlation_analysis.png - 相关性热图\n")
  cat("2. HOXC4_expression_comparison.png - 表达比较图\n")
  
}, error = function(e) {
  cat("执行过程中出现错误:", e$message, "\n")
})

# 示例2：使用自定义参数运行分析
cat("\n====================================\n")
cat("示例2：使用自定义参数运行分析\n")

# 如果你已经有一个真实的Seurat对象，你可以这样使用：
# results <- plot_gene_specific_correlation(
#   seurat_object = your_real_seurat_object,  # 替换为你的Seurat对象
#   gene_name = "你的目标基因",  # 例如 "DOCK5", "HOXC4"等
#   output_dir = "./your_output_directory",  # 自定义输出目录
#   top_percentile = 15,  # 自定义高表达阈值
#   bottom_percentile = 15,  # 自定义低表达阈值
#   assay = "RNA",  # 使用的assay
#   slot = "data",  # 使用的数据slot
#   plot_width = 14,  # 图片宽度
#   plot_height = 12,  # 图片高度
#   plot_dpi = 300  # 图片DPI
# )

cat("\n使用说明:\n")
cat("1. 准备好你的Seurat对象，确保其中包含目的基因的表达数据\n")
cat("2. 确保Seurat对象中包含细胞类型信息（列名为cell_type、seurat_clusters等）\n")
cat("3. 调用plot_gene_specific_correlation函数，传入必要的参数\n")
cat("4. 函数会自动根据基因表达水平分组，计算细胞通讯矩阵，生成相关性热图\n")
cat("5. 结果将保存在指定的输出目录中\n")
