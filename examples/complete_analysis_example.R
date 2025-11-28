# 完整分析示例脚本
# 展示如何使用 CellChatCompare 包进行基因表达分组和细胞通讯比较分析

# 加载必要的包
library(Seurat)
library(CellChat)
library(CellChatCompare)

# 设置工作目录
setwd("d:\\桌面\\celloracle\\m1\\CellChatCompare")

# 示例：创建模拟数据（实际使用时替换为真实数据）
# 这里我们创建一个简单的 Seurat 对象用于演示
create_simulated_data <- function() {
  # 创建模拟的基因表达矩阵
  n_cells <- 200
  n_genes <- 1000
  
  # 随机生成表达数据
  expr_matrix <- matrix(rpois(n_cells * n_genes, lambda = 0.5), nrow = n_genes, ncol = n_cells)
  
  # 添加一些差异表达的基因
  rownames(expr_matrix) <- paste0("Gene_", 1:n_genes)
  colnames(expr_matrix) <- paste0("Cell_", 1:n_cells)
  
  # 添加目标基因 DOCK5，其中50%的细胞高表达
  expr_matrix["DOCK5", 1:100] <- rpois(100, lambda = 5)  # 高表达组
  expr_matrix["DOCK5", 101:200] <- 0  # 低表达组
  
  # 创建细胞类型信息
  cell_types <- c(rep("CellType1", 50), rep("CellType2", 50), 
                  rep("CellType3", 50), rep("CellType4", 50))
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = "SimulatedData")
  seurat_obj$cell_subtype <- cell_types
  
  # 标准化和降维
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
  
  return(seurat_obj)
}

# 创建模拟数据
seurat_obj <- create_simulated_data()

# 查看数据基本信息
print(seurat_obj)
print(table(seurat_obj$cell_subtype))

# 1. 调用主函数进行完整分析
# 这里使用 auto 模式，当中位数为0时自动使用二进制分组
results <- compare_cellchat(
  seurat_obj = seurat_obj,          # Seurat 对象
  target_gene = "DOCK5",            # 目标基因
  celltype_col = "cell_subtype",    # 细胞类型列名
  group_method = "auto",            # 自动分组模式
  species = "human",                # 物种
  min_cells_high = 3,               # 高表达组最小细胞数
  min_cells_low = 10,               # 低表达组最小细胞数
  comm_level = "interaction",       # 通讯矩阵级别
  output_dir = "results_example",   # 结果输出目录
  plot_type = "all",                # 生成所有类型的图
  save_objects = TRUE               # 保存分析对象
)

# 2. 查看分析结果
cat("\n=== 分析结果摘要 ===\n")
cat("高表达组细胞数:", ncol(results$seurat_groups$high_expr), "\n")
cat("低表达组细胞数:", ncol(results$seurat_groups$low_expr), "\n")
cat("分组阈值:", results$seurat_groups$threshold, "\n")
cat("Mantel 测试结果 - 统计量:", results$mantel_result$statistic, "\n")
cat("Mantel 测试结果 - P值:", results$mantel_result$p_value, "\n")

# 3. 查看生成的图
cat("\n=== 生成的图 ===\n")
print(names(results$plots))

# 4. 单独查看某个图
# 如果需要在 RStudio 中查看图，可以使用以下命令
# print(results$plots$correlation)

# 5. 保存所有图为 PDF
pdf("results_example/all_plots.pdf", width = 12, height = 10)
for (plot_name in names(results$plots)) {
  print(results$plots[[plot_name]])
  title(main = plot_name)
}
dev.off()

cat("\n=== 分析完成 ===\n")
cat("结果已保存到:", "results_example", "\n")
