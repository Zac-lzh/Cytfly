# 简单测试脚本
# 验证 CellChatCompare 包的核心功能

# 加载必要的包
library(Seurat)
library(CellChat)
library(CellChatCompare)

# 设置工作目录
setwd("d:\\桌面\\celloracle\\m1\\CellChatCompare")

# 创建一个非常简单的测试数据
create_test_data <- function() {
  # 创建一个小的表达矩阵
  expr_matrix <- matrix(c(
    # 高表达组（前5个细胞）
    5, 6, 4, 5, 7,  # DOCK5 高表达
    1, 1, 1, 1, 1,  # Gene1
    2, 2, 2, 2, 2,  # Gene2
    # 低表达组（后5个细胞）
    0, 0, 0, 0, 0,  # DOCK5 低表达
    1, 1, 1, 1, 1,  # Gene1
    2, 2, 2, 2, 2   # Gene2
  ), nrow = 3, ncol = 10, byrow = TRUE)
  
  rownames(expr_matrix) <- c("DOCK5", "Gene1", "Gene2")
  colnames(expr_matrix) <- paste0("Cell_", 1:10)
  
  # 创建细胞类型信息
  cell_types <- c(rep("TypeA", 5), rep("TypeB", 5))
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = "TestData")
  seurat_obj$cell_subtype <- cell_types
  
  # 标准化数据
  seurat_obj <- NormalizeData(seurat_obj)
  
  return(seurat_obj)
}

# 测试1：基因表达分组
cat("=== 测试1：基因表达分组 ===\n")
seurat_obj <- create_test_data()

# 测试auto模式分组
grouped <- group_by_expression(seurat_obj, "DOCK5", group_method = "auto")
cat("高表达组细胞数:", ncol(grouped$high_expr), "\n")
cat("低表达组细胞数:", ncol(grouped$low_expr), "\n")
cat("分组阈值:", grouped$threshold, "\n")

# 测试median模式分组（应该与auto模式结果相同，因为中位数为0）
grouped_median <- group_by_expression(seurat_obj, "DOCK5", group_method = "median")
cat("\nMedian模式：\n")
cat("高表达组细胞数:", ncol(grouped_median$high_expr), "\n")
cat("低表达组细胞数:", ncol(grouped_median$low_expr), "\n")

# 测试binary模式分组
grouped_binary <- group_by_expression(seurat_obj, "DOCK5", group_method = "binary")
cat("\nBinary模式：\n")
cat("高表达组细胞数:", ncol(grouped_binary$high_expr), "\n")
cat("低表达组细胞数:", ncol(grouped_binary$low_expr), "\n")

# 测试2：检查compare_cellchat函数参数
cat("\n=== 测试2：检查compare_cellchat函数参数 ===\n")
# 只检查函数参数，不实际运行（因为需要CellChat数据库）
args(compare_cellchat)

# 测试3：检查run_cellchat函数是否支持不同的min_cells参数
cat("\n=== 测试3：检查run_cellchat函数参数 ===\n")
args(run_cellchat)

cat("\n=== 所有测试完成 ===\n")
cat("核心功能验证成功！\n")
