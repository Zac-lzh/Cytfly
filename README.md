# CellChatCompare: Comparative Analysis of Cell Communication

## Overview

CellChatCompare is an R package that provides tools to compare cell-cell communication networks between high and low expression groups of a target gene in single-cell RNA-seq data. It integrates Seurat, CellChat, and statistical analysis to identify differences in communication patterns.

## Installation

```r
# Install the package
devtools::install("path/to/CellChatCompare")

# Load the package
library(CellChatCompare)
```

## Quick Start

```r
# Load your Seurat object
seurat_obj <- readRDS("path/to/your/seurat_obj.rds")

# Run the complete analysis
result <- compare_cellchat(
  seurat_obj = seurat_obj,
  target_gene = "Hoxc4",
  celltype_col = "cell_type",
  group_method = "median",
  species = "mouse",
  output_dir = "cellchat_results"
)

# Visualize the results
plot_comparison(
  cellchat_high = result$cellchat_high,
  cellchat_low = result$cellchat_low,
  comparison_result = result$mantel_result,
  output_dir = "cellchat_results/plots"
)
```

## Detailed Usage

### 1. Group cells by target gene expression

```r
# Group cells by median expression
grouped_obj <- group_by_expression(
  seurat_obj = seurat_obj,
  target_gene = "Hoxc4",
  group_method = "median"
)

# Or group by percentile
grouped_obj <- group_by_expression(
  seurat_obj = seurat_obj,
  target_gene = "Hoxc4",
  group_method = "percentile",
  percentile = 0.75  # Top 25% as high expression
)
```

### 2. Run CellChat analysis

```r
# Run CellChat on high expression group
cellchat_high <- run_cellchat(
  seurat_obj = grouped_obj$high_expr,
  celltype_col = "cell_type",
  species = "mouse"
)

# Run CellChat on low expression group
cellchat_low <- run_cellchat(
  seurat_obj = grouped_obj$low_expr,
  celltype_col = "cell_type",
  species = "mouse"
)
```

### 3. Extract communication matrices

```r
# Extract pathway-level communication matrices
comm_high <- extract_communication_matrix(cellchat_high, "pathway")
comm_low <- extract_communication_matrix(cellchat_low, "pathway")

# Or extract interaction-level matrices
comm_high_int <- extract_communication_matrix(cellchat_high, "interaction")
comm_low_int <- extract_communication_matrix(cellchat_low, "interaction")
```

### 4. Perform Mantel test

```r
# Compare communication matrices using Mantel test
mantel_result <- mantel_test(
  matrix1 = comm_high,
  matrix2 = comm_low,
  method = "pearson",
  permutations = 999
)
```

### 5. Create visualizations

```r
# Create all plot types
plots <- plot_comparison(
  cellchat_high = cellchat_high,
  cellchat_low = cellchat_low,
  comparison_result = mantel_result,
  output_dir = "plots",
  plot_type = "all"
)

# Or create specific plot types
plots_hierarchy <- plot_comparison(
  cellchat_high = cellchat_high,
  cellchat_low = cellchat_low,
  output_dir = "plots",
  plot_type = "hierarchy"
)
```

## Parameters

### compare_cellchat()

| Parameter | Description | Default |
|-----------|-------------|---------|
| seurat_obj | A Seurat object | Required |
| target_gene | Target gene name | Required |
| celltype_col | Column name with cell type info | Required |
| group_method | Method to split cells ("median" or "percentile") | "median" |
| percentile | Percentile threshold | 0.5 |
| species | Species ("human" or "mouse") | "mouse" |
| assay | Assay to use | "RNA" |
| slot | Slot to use | "data" |
| min_cells | Min cells per cell type for CellChat | 10 |
| max_expr | Max expression threshold for CellChat | 0.95 |
| min_pct | Min percentage of cells expressing a gene | 0.1 |
| comm_level | Level of communication matrix ("pathway" or "interaction") | "pathway" |
| mantel_method | Correlation method for Mantel test | "pearson" |
| permutations | Number of permutations for significance testing | 999 |
| plot_type | Type of plot to generate ("all", "hierarchy", "circle", "heatmap") | "all" |
| max_links | Maximum number of links to display | 50 |
| output_dir | Directory to save results | "cellchat_results" |
| save_objects | Whether to save CellChat objects | TRUE |

## Output

The analysis generates the following outputs:

1. **Grouped Seurat objects**: High and low expression groups
2. **CellChat objects**: Communication analysis for each group
3. **Communication matrices**: Probability matrices for each group
4. **Mantel test results**: Statistical comparison of matrices
5. **Visualizations**: Various plots comparing communication networks
6. **Complete results object**: An R list containing all results

## Troubleshooting

### Common Issues

1. **Target gene not found**: Ensure the gene name exactly matches the rownames of your Seurat object
2. **Cell type column not found**: Check the column name in your Seurat metadata
3. **Insufficient cells**: Increase the `min_cells` parameter or ensure you have enough cells in each group
4. **Memory issues**: For large datasets, consider reducing the number of cells or using a subset

### Tips

- For large datasets, consider using a subset of cells for initial testing
- Ensure your Seurat object has normalized expression data
- Check that cell type annotations are accurate before running CellChat
- For human data, make sure to set `species = "human"`

## Citation

If you use this package, please cite:

- Seurat: Stuart et al., Cell, 2019
- CellChat: Jin et al., Nature Communications, 2021
- vegan: Oksanen et al., Community Ecology Package, 2022