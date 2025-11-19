#' @title Perform Mantel test between communication matrices
#' @description This function performs a Mantel test to compare two communication probability matrices
#' @param matrix1 First communication matrix
#' @param matrix2 Second communication matrix
#' @param method Correlation method ("pearson", "spearman", or "kendall")
#' @param permutations Number of permutations for significance testing
#' @param na.rm Logical indicating whether NA values should be removed
#' @return A list containing Mantel test statistics and p-value
#' @examples
#' mantel_result <- mantel_test(comm_matrix1, comm_matrix2, "pearson", 999)
#' @export
mantel_test <- function(matrix1, matrix2, method = "pearson", permutations = 999, na.rm = FALSE) {
  
  # Load vegan package for Mantel test
  library(vegan)
  
  # Check if matrices are square and have the same dimensions
  if (nrow(matrix1) != ncol(matrix1) || nrow(matrix2) != ncol(matrix2)) {
    stop("Both matrices must be square")
  }
  
  if (nrow(matrix1) != nrow(matrix2)) {
    stop("Matrices must have the same dimensions")
  }
  
  # Convert matrices to distance matrices (1 - probability)
  # This is because Mantel test works on distance matrices
  dist1 <- as.dist(1 - matrix1)
  dist2 <- as.dist(1 - matrix2)
  
  # Perform Mantel test
  mantel_result <- vegan::mantel(dist1, dist2, method = method, 
                                 permutations = permutations, na.rm = na.rm)
  
  # Extract results
  statistic <- mantel_result$statistic
  p_value <- mantel_result$signif
  
  # Create result object
  result <- list(
    statistic = statistic,
    p_value = p_value,
    method = method,
    permutations = permutations,
    na.rm = na.rm,
    mantel_obj = mantel_result
  )
  
  # Print results
  cat("Mantel test results:\n")
  cat("Method:", method, "\n")
  cat("Statistic:", round(statistic, 4), "\n")
  cat("P-value:", p_value, "\n")
  cat("Permutations:", permutations, "\n")
  
  if (p_value < 0.05) {
    cat("The communication matrices are significantly correlated (p < 0.05)\n")
  } else {
    cat("The communication matrices are not significantly correlated (p >= 0.05)\n")
  }
  
  return(result)
}