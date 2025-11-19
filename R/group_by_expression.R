#' @title Group cells by target gene expression
#' @description This function splits a Seurat object into high and low expression groups based on a target gene
#' @param seurat_obj A Seurat object
#' @param target_gene Character string specifying the target gene name
#' @param group_method Method to split cells ("median" or "percentile")
#' @param percentile Percentile threshold (only used when group_method = "percentile")
#' @param assay Assay to use (default "RNA")
#' @param slot Slot to use (default "data")
#' @return A list containing two Seurat objects: high_expr and low_expr
#' @examples
#' grouped_obj <- group_by_expression(seurat_obj, "Hoxc4", "median")
#' @export
group_by_expression <- function(seurat_obj, target_gene, group_method = "median", 
                                percentile = 0.5, assay = "RNA", slot = "data") {
  
  # Check if target gene exists in the object
  if (!target_gene %in% rownames(seurat_obj)) {
    stop(paste("Target gene", target_gene, "not found in the Seurat object"))
  }
  
  # Extract expression values
  expr_values <- FetchData(seurat_obj, vars = target_gene, assay = assay, slot = slot)
  expr_values <- expr_values[[target_gene]]
  
  # Determine threshold
  if (group_method == "median") {
    threshold <- median(expr_values)
  } else if (group_method == "percentile") {
    threshold <- quantile(expr_values, percentile)
  } else {
    stop("group_method must be either 'median' or 'percentile'")
  }
  
  # Create group labels
  group_labels <- ifelse(expr_values > threshold, "high", "low")
  
  # Add group information to metadata
  seurat_obj$expr_group <- factor(group_labels, levels = c("low", "high"))
  
  # Split the object
  high_expr <- subset(seurat_obj, subset = expr_group == "high")
  low_expr <- subset(seurat_obj, subset = expr_group == "low")
  
  # Print summary
  cat("Grouping cells by", target_gene, "expression using", group_method, "method\n")
  cat("Threshold:", threshold, "\n")
  cat("High expression group:", ncol(high_expr), "cells\n")
  cat("Low expression group:", ncol(low_expr), "cells\n")
  
  return(list(high_expr = high_expr, low_expr = low_expr, threshold = threshold))
}