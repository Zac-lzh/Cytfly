#' @title Run CellChat analysis on Seurat object
#' @description This function runs CellChat analysis on a Seurat object with improved error handling,
#'              progress tracking, and comprehensive validation
#' @param seurat_obj A Seurat object
#' @param celltype_col Column name in meta.data containing cell type information
#' @param species Species for CellChatDB ("human" or "mouse")
#' @param database Database to use ("CellChatDB" or custom)
#' @param min_cells Minimum number of cells per cell type (default 10)
#' @param max_expr Maximum expression threshold (default 0.95)
#' @param min_pct Minimum percentage of cells expressing a gene (default 0.1)
#' @param assay Assay name to use from Seurat object (default: "RNA")
#' @param slot Which slot to use for expression data (default: "data")
#' @param compute_pathways Whether to compute pathway-level communication (default: TRUE)
#' @param aggregate_network Whether to aggregate the network (default: TRUE)
#' @param ppi_projection Whether to perform PPI projection (default: TRUE)
#' @param verbose Whether to print detailed progress (default: TRUE)
#' @return A processed CellChat object with additional metadata and validation flags
#' @examples
#' cellchat_obj <- run_cellchat(seurat_obj, "cell_type", "mouse")
#' @export
run_cellchat <- function(seurat_obj, celltype_col, species = "mouse", 
                         database = "CellChatDB", min_cells = 10, 
                         max_expr = 0.95, min_pct = 0.1, assay = "RNA",
                         slot = "data", compute_pathways = TRUE,
                         aggregate_network = TRUE, ppi_projection = TRUE,
                         verbose = TRUE) {
  
  # 安全检查Seurat对象
  if (is.null(seurat_obj)) {
    stop("Error: Empty Seurat object provided")
  }
  
  # Check if required packages are installed
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("CellChat package is required but not installed. Please install it with devtools::install_github('sqjin/CellChat')")
  }
  
  # Load CellChat
  library(CellChat)
  
  # 安全检查细胞类型列
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", celltype_col, "not found in Seurat object metadata. Available columns:", 
               paste(colnames(seurat_obj@meta.data), collapse=", ")))
  }
  
  # Check if assay exists
  if (!(assay %in% names(seurat_obj@assays))) {
    stop(paste("Assay", assay, "not found in Seurat object. Available assays:", 
               paste(names(seurat_obj@assays), collapse=", ")))
  }
  
  # 初始化cellchat变量
  cellchat <- NULL
  
  # 创建进度追踪器
  progress <- list(
    steps = c("验证输入", "创建CellChat对象", "设置数据库", "数据子集", 
             "基因识别", "PPI投影", "计算通信概率", "过滤交互", 
             "通路水平计算", "网络聚合", "最终验证"),
    completed = 0,
    start_time = Sys.time()
  )
  
  # 进度输出函数
  print_progress <- function(step, success = TRUE) {
    if (verbose) {
      status <- ifelse(success, "✓", "✗")
      cat(status, step, "\n")
    }
  }
  
  # 验证输入完成
  print_progress("验证输入")
  
  # 创建CellChat对象与增强的错误处理
  if (verbose) cat("Creating CellChat object...\n")
  tryCatch({
    # 安全设置默认assay
    tryCatch({
      DefaultAssay(seurat_obj) <- assay
    }, error = function(e) {
      cat("Warning: Could not set default assay, continuing with current assay\n")
    })
    
    # 安全获取数据，支持新的slot参数
    data.input <- NULL
    meta <- NULL
    tryCatch({
      data.input <- GetAssayData(seurat_obj, slot = slot, assay = assay)
      meta <- seurat_obj@meta.data
      if (verbose) cat("Data extracted successfully - cells:", ncol(data.input), ", features:", nrow(data.input), "\n")
    }, error = function(e) {
      if (verbose) cat("Error extracting data:", e$message, "\n")
      # Try alternative way to get data
      tryCatch({
        if (isS4(seurat_obj) && "assays" %in% slotNames(seurat_obj) && 
            assay %in% names(seurat_obj@assays) && 
            slot %in% slotNames(seurat_obj@assays[[assay]])) {
          data.input <- as.matrix(seurat_obj@assays[[assay]]@data)
        } else {
          # 尝试通用方式获取数据
          stop("Cannot access data slot in Seurat object")
        }
        meta <- seurat_obj@meta.data
        if (verbose) cat("Data extracted using alternative method\n")
      }, error = function(e2) {
        stop("Failed to extract data from Seurat object: ", e2$message)
      })
    })
    
    # 创建CellChat对象
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = celltype_col)
    cat("✓ CellChat object created successfully\n")
  }, error = function(e) {
    cat("✗ Error creating CellChat object:", e$message, "\n")
    stop("Failed to create CellChat object: ", e$message)
  })
  
  # Set database
  cat("Setting CellChatDB for", species, "...\n")
  tryCatch({
    if (species == "human") {
      CellChatDB <- CellChat::CellChatDB.human
    } else if (species == "mouse") {
      CellChatDB <- CellChat::CellChatDB.mouse
    } else {
      stop("Species must be either 'human' or 'mouse'")
    }
    
    # Use a subset of CellChatDB for cell-cell communication analysis
    cellchat@DB <- CellChatDB
    cat("✓ Database set successfully\n")
  }, error = function(e) {
    cat("✗ Error setting CellChatDB:", e$message, "\n")
    # Fallback approach
    cat("Attempting to load database directly...\n")
    tryCatch({
      if (species == "human") {
        data(CellChatDB.human, package = "CellChat")
        CellChatDB <- CellChatDB.human
      } else {
        data(CellChatDB.mouse, package = "CellChat")
        CellChatDB <- CellChatDB.mouse
      }
      cellchat@DB <- CellChatDB
      cat("✓ Database loaded using fallback method\n")
    }, error = function(e2) {
      stop("Failed to load database: ", e2$message)
    })
  })
  
  # Subset the expression data
  cat("Subsetting data...\n")
  tryCatch({
    cellchat <- subsetData(cellchat)
    cat("✓ Data subset completed\n")
  }, error = function(e) {
    warning("Error during data subsetting: ", e$message)
  })
  
  # Identify over-expressed genes
  cat("Identifying over-expressed genes and interactions...\n")
  tryCatch({
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cat("✓ Gene identification completed\n")
  }, error = function(e) {
    warning("Error identifying over-expressed genes: ", e$message)
  })
  
  # Project gene expression data onto PPI network (可配置)
  if (ppi_projection) {
    if (verbose) cat("Projecting onto PPI network...\n")
    tryCatch({
      if (species == "human") {
        data(PPI.human, package = "CellChat")
        cellchat <- projectData(cellchat, PPI.human)
      } else {
        data(PPI.mouse, package = "CellChat")
        cellchat <- projectData(cellchat, PPI.mouse)
      }
      print_progress("PPI投影")
    }, error = function(e) {
      if (verbose) cat("✗ Warning: PPI projection failed, continuing without it:", e$message, "\n")
      print_progress("PPI投影", FALSE)
    })
  } else {
    if (verbose) cat("PPI projection skipped as requested\n")
    print_progress("PPI投影")
  }
  
  # Compute communication probability
  cat("Computing communication probability...\n")
  tryCatch({
    cellchat <- computeCommunProb(cellchat)
    cat("✓ Communication probability computed successfully\n")
  }, error = function(e) {
    cat("✗ Error computing communication probability:", e$message, "\n")
    # Try with different parameters
    cat("Trying with different parameters...\n")
    tryCatch({
      cellchat <- computeCommunProb(cellchat, type = "triMean")
      cat("✓ Communication probability computed with alternative parameters\n")
    }, error = function(e2) {
      stop("Failed to compute communication probability: ", e2$message)
    })
  })
  
  # Filter out insignificant interactions
  cat("Filtering interactions...\n")
  tryCatch({
    cellchat <- filterCommunication(cellchat, min.cells = min_cells)
    cat("✓ Interaction filtering completed\n")
  }, error = function(e) {
    warning("Error filtering interactions: ", e$message)
  })
  
  # Compute communication probability at pathway level (可配置)
  if (compute_pathways) {
    if (verbose) cat("Computing communication probability at pathway level...\n")
    tryCatch({
      cellchat <- computeCommunProbPathway(cellchat)
      print_progress("通路水平计算")
    }, error = function(e) {
      if (verbose) cat("✗ Warning: Pathway-level computation failed:", e$message, "\n")
      print_progress("通路水平计算", FALSE)
    })
  } else {
    if (verbose) cat("Pathway-level computation skipped as requested\n")
    print_progress("通路水平计算")
  }

  # Calculate aggregate network (可配置)
  if (aggregate_network) {
    if (verbose) cat("Aggregating communication network...\n")
    tryCatch({
      cellchat <- aggregateNet(cellchat)
      print_progress("网络聚合")
    }, error = function(e) {
      if (verbose) warning("Error aggregating network: ", e$message)
      print_progress("网络聚合", FALSE)
    })
  } else {
    if (verbose) cat("Network aggregation skipped as requested\n")
    print_progress("网络聚合")
  }
  
  # 辅助函数安全获取net字段
  get_net_safely <- function(obj) {
    if (isS4(obj) && "net" %in% slotNames(obj)) {
      return(obj@net)
    } else if (is.list(obj) && "net" %in% names(obj)) {
      return(obj$net)
    }
    return(NULL)
  }
  
  # 最终检查确保net字段存在且正确
  net <- get_net_safely(cellchat)
  if (is.null(net) || !"weight" %in% names(net) || is.null(net$weight)) {
    cat("CRITICAL: CellChat object does not have properly computed communication networks!\n")
    warning("CellChat analysis completed but with potential issues. Please check the results.")
    
    # 检查对象结构以帮助诊断
    if (isS4(cellchat)) {
      cat("CellChat object slots:", paste(slotNames(cellchat), collapse=", "), "\n")
    } else if (is.list(cellchat)) {
      cat("CellChat object fields:", paste(names(cellchat), collapse=", "), "\n")
    }
  } else {
    cat("✓ Communication networks computed successfully\n")
    weight_obj <- net$weight
    if (!is.null(weight_obj)) {
      if (is.matrix(weight_obj)) {
        cat("Communication network dimensions:", nrow(weight_obj), "x", ncol(weight_obj), "\n")
      } else if (is.list(weight_obj)) {
        cat("Communication network contains", length(weight_obj), "interactions\n")
      }
    }
  }
  
  # 最终验证完成
  print_progress("最终验证")
  
  # 计算执行时间
  execution_time <- difftime(Sys.time(), progress$start_time, units = "mins")
  
  # 添加分析元数据
  analysis_info <- list(
    n_cells = ncol(seurat_obj),
    n_celltypes = length(unique(seurat_obj@meta.data[[celltype_col]])),
    execution_time = execution_time,
    parameters = list(
      species = species,
      min_cells = min_cells,
      max_expr = max_expr,
      min_pct = min_pct,
      assay = assay,
      slot = slot,
      compute_pathways = compute_pathways,
      aggregate_network = aggregate_network,
      ppi_projection = ppi_projection
    ),
    success = !is.null(net) && !is.null(net$weight)
  )
  
  # 将元数据添加到结果对象
  if (isS4(cellchat)) {
    # 对于S4对象，尝试安全地添加元数据
    tryCatch({
      cellchat@analysis_info <- analysis_info
    }, error = function(e) {
      # 如果无法直接添加，返回带元数据的列表
      result <- list(cellchat = cellchat, analysis_info = analysis_info)
      if (verbose) cat("\nCellChat analysis completed for", ncol(seurat_obj), "cells and", 
                    length(unique(seurat_obj@meta.data[[celltype_col]])), "cell types\n")
      if (verbose) cat("Execution time:", round(execution_time, 2), "minutes\n")
      return(result)
    })
  } else if (is.list(cellchat)) {
    cellchat$analysis_info <- analysis_info
  }
  
  # 输出分析摘要
  if (verbose) {
    cat("\n=== CellChat Analysis Summary ===\n")
    cat("Cells analyzed:", ncol(seurat_obj), "\n")
    cat("Cell types:", length(unique(seurat_obj@meta.data[[celltype_col]])), "\n")
    cat("Execution time:", round(execution_time, 2), "minutes\n")
    if (!is.null(net) && !is.null(net$weight)) {
      cat("✓ Communication networks successfully computed\n")
    } else {
      cat("⚠ Partial results - some steps may have failed\n")
    }
    cat("================================\n")
  }
  
  return(cellchat)
}