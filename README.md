# 🔬 CellChatCompare 项目

## 📚 项目介绍

CellChatCompare 是一个用于比较和分析细胞通讯网络的 R 语言工具包，特别关注特定基因表达对细胞间通讯的影响。该工具基于 CellChat 包开发，提供了更丰富的可视化和分析功能。

## ✨ 主要功能

- **🔍 基因特异性细胞通讯分析**：根据特定基因的表达水平分组，分析不同表达组之间的细胞通讯差异
- **📊 相关性热图可视化**：使用 qcorrplot 生成美观的细胞通讯相关性热图，支持贝塞尔曲线连接
- **🧪 Mantel 测试**：进行细胞通讯网络的相关性检验
- **🎨 多种可视化类型**：支持层次图、圆形图、热图等多种可视化方式

## 📁 目录结构

```
CellChatCompare/
├── R/               # R 函数源码
├── examples/        # 示例脚本
├── output/          # 输出目录
├── DESCRIPTION      # 包描述文件
├── README.md        # 英文说明文档
└── README_中文.md   # 中文说明文档
```

## 📦 安装依赖

使用前请确保安装以下 R 包：

```r
install.packages(c("ggplot2", "patchwork", "corrr", "qgraph", "dplyr", "tidyr", "Seurat"))
# CellChat 包的安装
remotes::install_github("sqjin/CellChat")
```

## 🚀 使用示例

### 🔬 基因特异性相关性分析

```r
# 加载必要的包
library(Seurat)
library(CellChat)
library(ggplot2)

# 加载 Seurat 对象和 CellChat 对象
seurat_object <- readRDS("path/to/seurat_object.rds")
cellchat_high <- readRDS("path/to/cellchat_high.rds")
cellchat_low <- readRDS("path/to/cellchat_low.rds")

# 使用基因特异性相关性分析函数
source("R/plot_gene_specific_correlation.R")

# 分析特定基因对细胞通讯的影响
results <- plot_gene_specific_correlation(
  seurat_object = seurat_object,
  cellchat_high = cellchat_high,
  cellchat_low = cellchat_low,
  target_gene = "HOXC4",
  output_dir = "output/hoxc4_correlation"
)
```

## 🎯 核心函数

### plot_gene_specific_correlation

该函数是本项目的核心功能，用于分析特定基因表达与细胞通讯网络的关系。

**参数说明**：
- `seurat_object`：Seurat 对象，包含细胞的基因表达数据
- `cellchat_high`：高表达组的 CellChat 对象
- `cellchat_low`：低表达组的 CellChat 对象
- `target_gene`：要分析的目标基因名称
- `output_dir`：输出结果的目录路径
- `correlation_method`：相关性计算方法，默认为 "pearson"
- `pvalue_threshold`：显著性阈值，默认为 0.05

**返回值**：
- 包含相关性热图、表达比较图和 Mantel 测试结果的列表

## 📋 输入输出示例

#### 🔄 输入：
```r
# 调用函数示例
results <- plot_gene_specific_correlation(
  seurat_object = pbmc,
  cellchat_high = cellchat_high,
  cellchat_low = cellchat_low,
  target_gene = "HOXC4",
  output_dir = "output/hoxc4_result"
)
```

#### 📤 输出：
```
# 函数将生成以下文件：
# 1. output/hoxc4_result/correlation_heatmap.png - 细胞通讯相关性热图
# 2. output/hoxc4_result/gene_expression_comparison.png - 基因表达比较图
# 3. output/hoxc4_result/mantel_test_result.txt - Mantel 测试结果
```

## ⚠️ 注意事项

1. 确保 Seurat 对象中包含目标基因的表达数据
2. 对于大型数据集，可能需要调整参数以获得更好的可视化效果
3. 本工具支持不同版本的 Seurat 和 CellChat，但某些功能可能需要根据包版本进行调整

## 📝 许可证

本项目采用 MIT 许可证。

## 📬 联系方式

如有问题或建议，请通过以下方式联系：

- 邮箱：your.email@example.com
- GitHub：https://github.com/Zac-lzh/CellChatCompare