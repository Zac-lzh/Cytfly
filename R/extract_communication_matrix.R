#' @title Extract communication probability matrices from CellChat objects
#' @description This function extracts communication probability matrices from CellChat objects for comparison
#' @param cellchat_obj A CellChat object
#' @param level Level of communication matrix ("pathway" or "interaction")
#' @return A matrix of communication probabilities
#' @examples
#' comm_matrix <- extract_communication_matrix(cellchat_obj, "pathway")
#' @export
extract_communication_matrix <- function(cellchat_obj, level = "pathway") {
  # 检查对象是否为空
  if (is.null(cellchat_obj)) {
    cat("Error: Empty CellChat object provided\n")
    return(NULL)
  }
  
  # 安全获取net字段 - 支持S4对象和列表
  get_net <- function(obj) {
    if (isS4(obj)) {
      # 对于S4对象使用@访问器
      if ("net" %in% slotNames(obj)) {
        return(obj@net)
      }
    } else if (is.list(obj)) {
      # 对于列表使用$访问器
      if ("net" %in% names(obj)) {
        return(obj$net)
      }
    }
    return(NULL)
  }
  
  net <- get_net(cellchat_obj)
  
  # Check if CellChat object has been processed, and automatically fix if needed
  if (is.null(net) || !"weight" %in% names(net) || is.null(net$weight)) {
    cat("CellChat object not fully processed. Automatically running computeCommunProb...\n")
    tryCatch({
      cellchat_obj <- computeCommunProb(cellchat_obj)
      cat("computeCommunProb completed.\n")
      net <- get_net(cellchat_obj)
    }, error = function(e) {
      cat("Warning: Failed to run computeCommunProb:", e$message, "\n")
    })
  }
  
  # For pathway level, check and compute if needed
  if (level == "pathway") {
    if (is.null(net) || !"netP" %in% names(net) || is.null(net$netP) || !"weight" %in% names(net$netP)) {
      cat("Pathway-level communication not found. Running computeCommunProbPathway...\n")
      tryCatch({
        cellchat_obj <- computeCommunProbPathway(cellchat_obj)
        cat("computeCommunProbPathway completed.\n")
        net <- get_net(cellchat_obj)
      }, error = function(e) {
        cat("Warning: Failed to run computeCommunProbPathway:", e$message, "\n")
      })
    }
  }
  
  # Extract matrix based on level with improved error handling
  comm_matrix <- NULL
  
  # First try pathway level if requested
  if (level == "pathway" && !is.null(net) && "netP" %in% names(net) && !is.null(net$netP)) {
    if ("weight" %in% names(net$netP) && !is.null(net$netP$weight)) {
      comm_matrix <- net$netP$weight
      cat("Successfully extracted pathway-level communication matrix\n")
    } else {
      warning("Pathway-level communication matrix not available after attempted computation")
      level <- "interaction"  # Fall back to interaction level
    }
  } else if (level == "pathway") {
    warning("Pathway-level communication matrix not available, falling back to interaction level")
    level <- "interaction"
  }
  
  # Always try interaction level as fallback
  if (level == "interaction" || is.null(comm_matrix)) {
    if (!is.null(net)) {
      # First try weight field (processed probabilities)
      if ("weight" %in% names(net) && !is.null(net$weight)) {
        comm_matrix <- net$weight
        cat("Successfully extracted interaction-level communication matrix (weight)\n")
      } else if ("prob" %in% names(net) && !is.null(net$prob)) {
        # Fall back to prob field (raw probabilities)
        cat("Warning: Using prob field instead of weight (raw probabilities)\n")
        comm_matrix <- net$prob
      } else {
        # Final check: look for other possible communication matrices
        cat("Error: Cannot find communication matrix in CellChat object\n")
        if (length(names(net)) > 0) {
          cat("Available fields in net:", paste(names(net), collapse=", "), "\n")
        } else {
          cat("Net field is empty or not properly initialized\n")
        }
        return(NULL)  # 安全返回NULL而不是抛出错误
      }
    } else {
      cat("Error: No net field found in CellChat object\n")
      return(NULL)
    }
  }
  
  # Ensure matrix is square (sender x receiver)
  if (nrow(comm_matrix) != ncol(comm_matrix)) {
    cat("Warning: Communication matrix is not square. Taking the square portion.\n")
    min_dim <- min(nrow(comm_matrix), ncol(comm_matrix))
    comm_matrix <- comm_matrix[1:min_dim, 1:min_dim]
  }
  
  cat("Extracted", level, "-level communication matrix with dimensions:", 
      nrow(comm_matrix), "x", ncol(comm_matrix), "\n")
  
  return(comm_matrix)
}