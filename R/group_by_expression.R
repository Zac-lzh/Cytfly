#' @title Group cells by target gene expression
#' @description This function splits a Seurat object into high and low expression groups based on a target gene
#' @param seurat_obj A Seurat object
#' @param target_gene Character string specifying the target gene name
#' @param group_method Method to split cells ("median", "percentile", or "auto")
#' @param percentile Percentile threshold (only used when group_method = "percentile")
#' @param assay Assay to use (default "RNA")
#' @param slot Slot to use (default "data")
#' @return A list containing two Seurat objects: high_expr and low_expr
#' @examples
#' grouped_obj <- group_by_expression(seurat_obj, "DOCK5", "auto")
#' @export
group_by_expression <- function(seurat_obj, target_gene, group_method = "auto", 
                                percentile = 0.5, assay = "RNA", slot = "data") {
  
  # Check if target gene exists in the object
  if (!target_gene %in% rownames(seurat_obj)) {
    stop(paste("Target gene", target_gene, "not found in the Seurat object"))
  }
  
  # Extract expression values
  expr_values <- FetchData(seurat_obj, vars = target_gene, assay = assay, slot = slot)
  expr_values <- expr_values[[target_gene]]
  
  # Determine threshold and groups
  if (group_method == "auto") {
    # Auto mode: use median, but if median is 0, use binary grouping
    threshold <- median(expr_values)
    if (threshold == 0) {
      # Binary grouping: high (> 0), low (== 0)
      group_labels <- ifelse(expr_values > 0, "high", "low")
      cat("Median is 0, using binary grouping: high (> 0), low (== 0)\n")
    } else {
      # Median grouping
      group_labels <- ifelse(expr_values > threshold, "high", "low")
      cat("Using median grouping with threshold:", threshold, "\n")
    }
  } else if (group_method == "median") {
    threshold <- median(expr_values)
    # Check if median is 0, if so use binary grouping
    if (threshold == 0) {
      group_labels <- ifelse(expr_values > 0, "high", "low")
      cat("Median is 0, using binary grouping: high (> 0), low (== 0)\n")
    } else {
      group_labels <- ifelse(expr_values > threshold, "high", "low")
    }
  } else if (group_method == "percentile") {
    threshold <- quantile(expr_values, percentile)
    group_labels <- ifelse(expr_values > threshold, "high", "low")
  } else if (group_method == "binary") {
    # Explicit binary grouping
    threshold <- 0
    group_labels <- ifelse(expr_values > 0, "high", "low")
  } else {
    stop("group_method must be either 'median', 'percentile', 'binary', or 'auto'")
  }
  
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