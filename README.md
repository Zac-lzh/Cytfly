# 🔬 CellChatCompare

## 📚 项目介绍

CellChatCompare 是一个用于比较和分析细胞通讯网络的 R 语言工具包，特别关注特定基因表达对细胞间通讯的影响。该工具基于 CellChat 包开发，提供了从基因表达分组到细胞通讯分析再到结果可视化的完整工作流。

## ✨ 主要功能

- **🔍 基因表达分组**：根据目标基因的表达水平将细胞分为高表达组和低表达组
- **📊 CellChat 分析**：对高低表达组分别进行细胞通讯网络分析
- **🧪 Mantel 测试**：比较两组细胞通讯网络的相关性
- **🎨 多种可视化**：生成层次图、圆形图、热图和相关性图等多种可视化结果
- **📈 通讯矩阵对比**：直接比较高低表达组的细胞通讯矩阵

## 📁 目录结构

```
CellChatCompare/
├── R/               # R 函数源码
│   ├── compare_cellchat.R              # 主函数，整合所有功能
│   ├── extract_communication_matrix.R  # 提取通讯矩阵
│   ├── group_by_expression.R           # 按基因表达分组
│   ├── mantel_test.R                   # Mantel 测试
│   ├── plot_comparison.R               # 结果可视化
│   ├── plot_gene_specific_correlation.R # 基因特异性相关性分析
│   └── run_cellchat.R                  # 运行 CellChat 分析
├── examples/        # 示例脚本
├── .gitignore       # Git 忽略文件
├── DESCRIPTION      # 包描述文件
└── README.md        # 中文说明文档
```

## 📦 安装依赖

使用前请确保安装以下 R 包：

```r
# 安装 Seurat
devtools::install_github('satijalab/seurat', ref = 'seurat5')

# 安装 CellChat
devtools::install_github('sqjin/CellChat')

# 安装其他依赖包
install.packages(c("dplyr", "ggplot2", "reshape2", "vegan", "igraph", "patchwork", "RColorBrewer", "magrittr", "corrr"))

# 安装 linkET (用于 qcorrplot 函数)
install.packages("linkET")
```

## 🚀 使用示例

### 🔬 完整工作流示例

```r
# 加载必要的包
library(Seurat)
library(CellChat)
library(CellChatCompare)

# 假设您已经有一个 Seurat 对象 (seurat_obj)
# seurat_obj <- readRDS("path/to/seurat_object.rds")

# 1. 调用主函数进行完整分析
results <- compare_cellchat_by_gene(
  seurat_obj = seurat_obj,          # Seurat 对象
  target_gene = "DOCK5",            # 目标基因
  species = "human",                # 物种（"human" 或 "mouse"）
  celltype_col = "cell_subtype",    # 细胞类型列名
  group_method = "median",          # 分组方法（"median" 或 "percentile"）
  min_cells_high = 3,               # 高表达组最小细胞数
  min_cells_low = 10,               # 低表达组最小细胞数
  output_dir = "results",           # 结果输出目录
  plot_type = "all"                 # 绘制所有类型的图
)

# 2. 查看结果
names(results)  # 查看生成的图
results$correlation  # 查看相关性图
```

### 📊 分步使用示例

```r
# 加载必要的包
library(Seurat)
library(CellChat)
library(CellChatCompare)

# 1. 按基因表达分组
grouped <- group_by_expression(
  seurat_obj = seurat_obj, 
  target_gene = "DOCK5", 
  group_method = "median"
)
seurat_high <- grouped$high_expr
seurat_low <- grouped$low_expr

# 2. 运行 CellChat 分析
cellchat_high <- run_cellchat(
  seurat_obj = seurat_high, 
  celltype_col = "cell_subtype", 
  species = "human"
)

cellchat_low <- run_cellchat(
  seurat_obj = seurat_low, 
  celltype_col = "cell_subtype", 
  species = "human"
)

# 3. 提取通讯矩阵
comm_high <- extract_communication_matrix(cellchat_high, "interaction")
comm_low <- extract_communication_matrix(cellchat_low, "interaction")

# 4. 进行 Mantel 测试
mantel_result <- mantel_test(comm_high, comm_low, method = "pearson")

# 5. 绘制比较图
plots <- plot_comparison(
  cellchat_high = cellchat_high, 
  cellchat_low = cellchat_low, 
  comparison_result = mantel_result,
  output_dir = "results",
  gene_name = "DOCK5"
)
```

## 🎯 核心函数

### compare_cellchat_by_gene

**功能**：主函数，整合基因表达分组、CellChat 分析、Mantel 测试和结果可视化的完整工作流

**参数**：
- `seurat_obj`：Seurat 对象，包含细胞的基因表达数据
- `target_gene`：目标基因名称
- `species`：物种（"human" 或 "mouse"）
- `celltype_col`：细胞类型列名
- `group_method`：分组方法（"median" 或 "percentile"）
- `percentile`：百分位数阈值（仅当 `group_method = "percentile"` 时使用）
- `min_cells_high`：高表达组最小细胞数
- `min_cells_low`：低表达组最小细胞数
- `output_dir`：结果输出目录
- `plot_type`：绘制的图类型（"all"、"hierarchy"、"circle"、"heatmap" 或 "correlation"）
- `max_links`：最大连接数
- `with_couple`：是否添加连接线

**返回值**：包含所有分析结果和图的列表

### group_by_expression

**功能**：根据目标基因的表达水平将细胞分为高表达组和低表达组

**参数**：
- `seurat_obj`：Seurat 对象
- `target_gene`：目标基因名称
- `group_method`：分组方法（"median" 或 "percentile"）
- `percentile`：百分位数阈值
- `assay`：使用的 assay
- `slot`：使用的数据槽

**返回值**：包含高表达组、低表达组和阈值的列表

### run_cellchat

**功能**：对 Seurat 对象运行 CellChat 分析

**参数**：
- `seurat_obj`：Seurat 对象
- `celltype_col`：细胞类型列名
- `species`：物种
- `database`：使用的数据库
- `min_cells`：最小细胞数
- `max_expr`：最大表达阈值
- `min_pct`：最小表达百分比
- `assay`：使用的 assay
- `slot`：使用的数据槽
- `compute_pathways`：是否计算通路水平通讯
- `aggregate_network`：是否聚合网络
- `ppi_projection`：是否进行 PPI 投影
- `verbose`：是否显示详细进度

**返回值**：处理后的 CellChat 对象

### mantel_test

**功能**：比较两个通讯矩阵的相关性

**参数**：
- `matrix1`：第一个通讯矩阵
- `matrix2`：第二个通讯矩阵
- `method`：相关性方法（"pearson"、"spearman" 或 "kendall"）
- `permutations`：置换次数
- `na.rm`：是否移除 NA 值

**返回值**：包含 Mantel 测试结果的列表

### plot_comparison

**功能**：绘制通讯网络比较图

**参数**：
- `cellchat_high`：高表达组的 CellChat 对象
- `cellchat_low`：低表达组的 CellChat 对象
- `comparison_result`：Mantel 测试结果
- `output_dir`：输出目录
- `plot_type`：图类型
- `max_links`：最大连接数
- `with_couple`：是否添加连接线
- `gene_name`：基因名称

**返回值**：包含所有图的列表

## 📋 输入输出示例

#### 🔄 输入：
```r
# 调用主函数示例
results <- compare_cellchat_by_gene(
  seurat_obj = seurat_obj,
  target_gene = "DOCK5",
  species = "human",
  celltype_col = "cell_subtype",
  output_dir = "results"
)
```

#### 📤 输出：
```
# 函数将生成以下文件：
# 1. results/communication_hierarchy_comparison.png - 层次图比较
# 2. results/communication_circle_comparison.png - 圆形图比较
# 3. results/communication_difference_heatmap.png - 差异热图
# 4. results/communication_correlation_comparison.png - 相关性图
# 5. results/mantel_test_summary.png - Mantel 测试结果
```

## ⚠️ 注意事项

1. 确保 Seurat 对象中包含目标基因的表达数据
2. 确保细胞类型列名正确，且包含有效的细胞类型信息
3. 对于大型数据集，可能需要调整 `min_cells` 参数以获得更好的结果
4. 建议使用 Seurat 5 和最新版本的 CellChat 以获得最佳兼容性
5. 运行时间可能较长，具体取决于数据集大小和计算资源

## 📝 许可证

本项目采用 MIT 许可证。

## 📬 联系方式

如有问题或建议，请通过以下方式联系：

- 邮箱：18960339395@163.com
- GitHub：https://github.com/Zac-lzh/CellChatCompare
