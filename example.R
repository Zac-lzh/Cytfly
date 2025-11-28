# Example usage of CellChatCompare package

# Load required packages
library(Seurat)
library(CellChatCompare)

# Load your Seurat object
seurat_obj <- readRDS("d:/桌面/celloracle/m1/seurat_obj.rds")

# Check the structure of your Seurat object
print(seurat_obj)
print(head(seurat_obj@meta.data))

# Run the complete analysis
result <- compare_cellchat(
  seurat_obj = seurat_obj,
  target_gene = "Hoxc4",
  celltype_col = "cellType",  # Replace with your actual cell type column name
  group_method = "median",
  species = "mouse",
  output_dir = "d:/桌面/celloracle/m1/Hoxc4_cellchat_results"
)

# Explore the results
print(names(result))

# Access the Mantel test results
print(result$mantel_result)

# View the plots
print(result$plots$hierarchy)
print(result$plots$circle)
print(result$plots$heatmap)

# If you want to run the analysis step by step:

# 1. Group cells by expression
grouped_obj <- group_by_expression(
  seurat_obj = seurat_obj,
  target_gene = "Hoxc4",
  group_method = "median"
)

# 2. Run CellChat on each group
cellchat_high <- run_cellchat(
  seurat_obj = grouped_obj$high_expr,
  celltype_col = "cell_type",
  species = "mouse"
)

cellchat_low <- run_cellchat(
  seurat_obj = grouped_obj$low_expr,
  celltype_col = "cell_type",
  species = "mouse"
)

# 3. Extract communication matrices
comm_high <- extract_communication_matrix(cellchat_high, "pathway")
comm_low <- extract_communication_matrix(cellchat_low, "pathway")

# 4. Perform Mantel test
mantel_result <- mantel_test(comm_high, comm_low, "pearson", 999)

# 5. Create visualizations
plots <- plot_comparison(
  cellchat_high = cellchat_high,
  cellchat_low = cellchat_low,
  comparison_result = mantel_result,
  output_dir = "d:/桌面/celloracle/m1/Hoxc4_plots",
  plot_type = "all"
)