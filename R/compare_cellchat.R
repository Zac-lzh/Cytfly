#' @title Compare cell-cell communication between high and low expression groups
#' @description This function performs the complete analysis pipeline to compare cell-cell 
#' communication between high and low expression groups of a target gene
#' @param seurat_obj A Seurat object
#' @param target_gene Character string specifying the target gene name
#' @param celltype_col Column name in metadata containing cell type information
#' @param group_method Method to split cells ("median", "percentile", "binary", or "auto")
#' @param percentile Percentile threshold (only used when group_method = "percentile")
#' @param species Species ("human" or "mouse")
#' @param assay Assay to use (default "RNA")
#' @param slot Slot to use (default "data")
#' @param max_expr Maximum expression threshold for CellChat (default 0.95)
#' @param min_pct Minimum percentage of cells expressing a gene for CellChat (default 0.1)
#' @param comm_level Level of communication matrix ("pathway" or "interaction")
#' @param mantel_method Correlation method for Mantel test ("pearson", "spearman", or "kendall")
#' @param permutations Number of permutations for Mantel test significance testing
#' @param plot_type Type of plot to generate ("all", "hierarchy", "circle", "heatmap", or "correlation")
#' @param max_links Maximum number of links to display in plots
#' @param save_objects Whether to save CellChat objects (default TRUE)
#' @param with_couple Logical, whether to add connection lines in correlation plot
#' @return A list containing all analysis results
#' @examples
#' result <- compare_cellchat(
#'   seurat_obj = seurat_obj,
#'   target_gene = "DOCK5",
#'   celltype_col = "cell_subtype",
#'   group_method = "auto",
#'   species = "human"
#' )
#' @export
compare_cellchat <- function(seurat_obj, target_gene, celltype_col, 
                             group_method = "auto", percentile = 0.5, 
                             species = "human", assay = "RNA", slot = "data",
                             max_expr = 0.95, min_pct = 0.1,
                             comm_level = "interaction", mantel_method = "pearson", 
                             permutations = 999, plot_type = "all", max_links = 50,
                             save_objects = TRUE,
                             with_couple = TRUE) {
  
  # 设置默认参数
  min_cells_high <- 3  # 高表达组最小细胞数
  min_cells_low <- 10  # 低表达组最小细胞数
  output_dir <- paste0("cellchat_results_", target_gene)  # 结果输出目录
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Starting CellChat comparison analysis for gene:", target_gene, "\n")
  cat("Output directory:", output_dir, "\n")
  
  # Step 1: Group cells by target gene expression
  cat("\n=== Step 1: Grouping cells by expression ===\n")
  grouped_obj <- group_by_expression(seurat_obj, target_gene, 
                                      group_method, percentile, assay, slot)
  
  # Save grouped objects if requested
  if (save_objects) {
    saveRDS(grouped_obj$high_expr, file.path(output_dir, "high_expression_seurat.rds"))
    saveRDS(grouped_obj$low_expr, file.path(output_dir, "low_expression_seurat.rds"))
    cat("Saved grouped Seurat objects\n")
  }
  
  # Step 2: Run CellChat analysis on each group
  cat("\n=== Step 2: Running CellChat analysis ===\n")
  cat("Analyzing high expression group...\n")
  cellchat_high <- run_cellchat(grouped_obj$high_expr, celltype_col, species, 
                                 min_cells = min_cells_high, max_expr = max_expr, min_pct = min_pct)
  
  cat("Analyzing low expression group...\n")
  cellchat_low <- run_cellchat(grouped_obj$low_expr, celltype_col, species, 
                                min_cells = min_cells_low, max_expr = max_expr, min_pct = min_pct)
  
  # Save CellChat objects if requested
  if (save_objects) {
    saveRDS(cellchat_high, file.path(output_dir, "high_expression_cellchat.rds"))
    saveRDS(cellchat_low, file.path(output_dir, "low_expression_cellchat.rds"))
    cat("Saved CellChat objects\n")
  }
  
  # Step 3: Extract communication matrices
  cat("\n=== Step 3: Extracting communication matrices ===\n")
  comm_high <- extract_communication_matrix(cellchat_high, comm_level)
  comm_low <- extract_communication_matrix(cellchat_low, comm_level)
  
  # Save communication matrices
  write.csv(comm_high, file.path(output_dir, "high_expression_comm_matrix.csv"))
  write.csv(comm_low, file.path(output_dir, "low_expression_comm_matrix.csv"))
  
  # Step 4: Perform Mantel test
  cat("\n=== Step 4: Performing Mantel test ===\n")
  mantel_result <- mantel_test(comm_high, comm_low, mantel_method, permutations)
  
  # Save Mantel test results
  write.csv(data.frame(
    Statistic = mantel_result$statistic,
    P_value = mantel_result$p_value,
    Method = mantel_result$method,
    Permutations = mantel_result$permutations
  ), file.path(output_dir, "mantel_test_results.csv"), row.names = FALSE)
  
  # Step 5: Create visualizations
  cat("\n=== Step 5: Creating visualizations ===\n")
  plots <- plot_comparison(cellchat_high, cellchat_low, mantel_result, 
                           output_dir, plot_type, max_links,
                           with_couple = with_couple, gene_name = target_gene)
  
  # Step 6: Compile results
  results <- list(
    seurat_groups = grouped_obj,
    cellchat_high = cellchat_high,
    cellchat_low = cellchat_low,
    comm_high = comm_high,
    comm_low = comm_low,
    mantel_result = mantel_result,
    plots = plots,
    parameters = list(
      target_gene = target_gene,
      celltype_col = celltype_col,
      group_method = group_method,
      percentile = percentile,
      species = species,
      min_cells_high = min_cells_high,
      min_cells_low = min_cells_low,
      comm_level = comm_level,
      mantel_method = mantel_method,
      permutations = permutations,
      output_dir = output_dir
    )
  )
  
  # Save results object
  if (save_objects) {
    saveRDS(results, file.path(output_dir, "complete_analysis_results.rds"))
    cat("Saved complete analysis results\n")
  }
  
  cat("\n=== Analysis completed successfully! ===\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(results)
}