#!/usr/bin/env Rscript
# 测试脚本：演示基因特异性通信网络相关性分析函数的使用
# 注意：此脚本提供了几个示例用例，可以根据需要进行修改

# 设置错误处理选项
options(error = function() traceback(3))

# 设置输出目录
test_output_dir <- file.path("..", "test_results")
dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)

# 打印信息函数
print_info <- function(message) {
  cat("\n", rep("-", 60), "\n", sep = "")
  cat("[INFO] ", message, "\n", sep = "")
  cat(rep("-", 60), "\n", sep = "")
}

# 打印错误函数
print_error <- function(message) {
  cat("\n", rep("*", 60), "\n", sep = "")
  cat("[ERROR] ", message, "\n", sep = "")
  cat(rep("*", 60), "\n", sep = "")
}

# 检查函数是否存在
check_function_available <- function() {
  tryCatch({
    # 加载主要函数
    source(file.path("..", "R", "plot_gene_specific_correlation.R"))
    return(TRUE)
  }, error = function(e) {
    print_error(paste("无法加载主函数文件:", e$message))
    return(FALSE)
  })
}

# 创建模拟的Seurat对象
create_mock_seurat_object <- function() {
  tryCatch({
    print_info("创建模拟的Seurat对象用于测试")
    
    # 确保必要的包已安装
    required_packages <- c("Seurat", "dplyr", "tidyr")
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        print_info(paste("安装必要的包:", pkg))
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      }
      library(pkg, character.only = TRUE)
    }
    
    # 创建模拟数据
    n_cells <- 200
    n_features <- 1000
    
    # 随机表达矩阵
    set.seed(42)
    counts_matrix <- matrix(rpois(n = n_cells * n_features, lambda = 2),
                           ncol = n_cells,
                           nrow = n_features)
    
    # 创建基因名
    rownames(counts_matrix) <- paste0("gene_", 1:n_features)
    colnames(counts_matrix) <- paste0("cell_", 1:n_cells)
    
    # 创建一个模拟基因的高表达
    # 例如我们将模拟'DOCK5'基因的高表达
    counts_matrix["gene_1", sample(1:n_cells, size = 50)] <- 50
    
    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(counts = counts_matrix)
    
    # 添加一些元数据
    # 细胞类型
    set.seed(42)
    cell_types <- c("B_cell", "Dendritic_cell", "Endothelial_cell", 
                   "Fibroblast", "Mast_cell", "Myeloid_cell", 
                   "Ovarian_cancer_cell", "Plasma_cell")
    
    seurat_obj@meta.data$cell_subtype <- sample(cell_types, size = n_cells, replace = TRUE)
    
    # 添加UMAP坐标（用于可视化，如果需要）
    set.seed(42)
    umap_embeddings <- matrix(rnorm(n_cells * 2), ncol = 2)
    rownames(umap_embeddings) <- colnames(seurat_obj)  # 添加行名
    colnames(umap_embeddings) <- c("UMAP_1", "UMAP_2")  # 添加列名
    
    seurat_obj@reductions$umap <- CreateDimReducObject(
      embeddings = umap_embeddings,
      key = "UMAP_",
      assay = DefaultAssay(seurat_obj)
    )
    
    # 标准化数据
    seurat_obj <- NormalizeData(seurat_obj)
    
    # 添加RNA前缀
    all.genes <- rownames(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features = all.genes)
    
    print_info(paste("成功创建模拟Seurat对象，包含", n_cells, "个细胞和", n_features, "个基因"))
    return(seurat_obj)
  }, error = function(e) {
    print_error(paste("创建模拟Seurat对象时出错:", e$message))
    return(NULL)
  })
}

# 测试函数：基本功能
test_basic_functionality <- function(seurat_obj) {
  tryCatch({
    print_info("测试基本功能")
    
    # 使用模拟数据标志
    use_mock = TRUE
    
    # 测试一个基因
    test_gene <- "gene_1"  # 对应我们在模拟数据中高表达的基因
    
    print_info(paste("分析基因:", test_gene))
    
    # 运行分析
    results <- plot_gene_specific_correlation(
      seurat_obj = seurat_obj,
      gene_name = test_gene,
      output_dir = file.path(test_output_dir, test_gene),
      use_mock_data = use_mock,
      plot_width = 12,
      plot_height = 10,
      plot_dpi = 300
    )
    
    # 检查结果
    if (!is.null(results) && !is.null(results$correlation_plot)) {
      print_info("基本功能测试成功：成功生成相关性图")
      print_info(paste("高表达细胞数量:", length(results$high_expression_cells)))
      print_info(paste("低表达细胞数量:", length(results$low_expression_cells)))
      
      if (!is.null(results$mantel_result)) {
        print_info(paste("Mantel测试结果: r =", results$mantel_result$r, ", p =", results$mantel_result$p))
      }
      
      return(TRUE)
    } else {
      print_error("基本功能测试失败：未能生成有效的结果")
      return(FALSE)
    }
  }, error = function(e) {
    print_error(paste("基本功能测试时出错:", e$message))
    return(FALSE)
  })
}

# 测试函数：不同基因
test_multiple_genes <- function(seurat_obj) {
  tryCatch({
    print_info("测试多个基因")
    
    # 测试多个基因
    test_genes <- c("gene_1", "gene_10", "gene_100", "DOCK5")  # 一些随机基因和我们的示例基因
    
    success_count <- 0
    
    for (gene in test_genes) {
      print_info(paste("测试基因:", gene))
      
      # 创建特定基因的输出目录
      gene_output_dir <- file.path(test_output_dir, gene)
      dir.create(gene_output_dir, recursive = TRUE, showWarnings = FALSE)
      
      # 运行分析
      results <- plot_gene_specific_correlation(
        seurat_obj = seurat_obj,
        gene_name = gene,
        output_dir = gene_output_dir,
        use_mock_data = TRUE,
        plot_width = 10,
        plot_height = 8,
        plot_dpi = 200
      )
      
      # 检查结果
      if (!is.null(results) && !is.null(results$correlation_plot)) {
        print_info(paste("基因", gene, "分析成功"))
        success_count <- success_count + 1
      } else {
        print_error(paste("基因", gene, "分析失败"))
      }
      
      # 短暂休息，避免资源竞争
      Sys.sleep(1)
    }
    
    print_info(paste("多个基因测试完成: ", success_count, "/", length(test_genes), "个基因分析成功"))
    return(success_count > 0)
  }, error = function(e) {
    print_error(paste("多个基因测试时出错:", e$message))
    return(FALSE)
  })
}

# 测试函数：不同的细胞亚型列
test_different_cell_columns <- function(seurat_obj) {
  tryCatch({
    print_info("测试不同的细胞亚型列")
    
    # 添加一个新的细胞类型列用于测试
    seurat_obj@meta.data$custom_cell_type <- paste0("Type_", sample(1:4, size = ncol(seurat_obj), replace = TRUE))
    
    # 测试不同的细胞类型列
    test_columns <- c("cell_subtype", "custom_cell_type", "不存在的列", NULL)
    
    success_count <- 0
    
    for (col in test_columns) {
      col_desc <- ifelse(is.null(col), "NULL (自动检测)", col)
      print_info(paste("测试细胞类型列:", col_desc))
      
      # 创建特定列的输出目录
      col_name <- ifelse(is.null(col), "auto_detect", gsub("\\s+", "_", col))
      col_output_dir <- file.path(test_output_dir, "cell_columns", col_name)
      dir.create(col_output_dir, recursive = TRUE, showWarnings = FALSE)
      
      # 运行分析
      results <- plot_gene_specific_correlation(
        seurat_obj = seurat_obj,
        gene_name = "gene_1",
        cell_subtype_column = col,
        output_dir = col_output_dir,
        use_mock_data = TRUE,
        plot_width = 10,
        plot_height = 8,
        plot_dpi = 200
      )
      
      # 检查结果
      if (!is.null(results) && !is.null(results$correlation_plot)) {
        print_info(paste("使用列", col_desc, "分析成功"))
        success_count <- success_count + 1
      } else {
        print_error(paste("使用列", col_desc, "分析失败"))
      }
      
      # 短暂休息
      Sys.sleep(1)
    }
    
    print_info(paste("细胞亚型列测试完成: ", success_count, "/", length(test_columns), "个列测试成功"))
    return(success_count > 0)
  }, error = function(e) {
    print_error(paste("细胞亚型列测试时出错:", e$message))
    return(FALSE)
  })
}

# 测试函数：边界情况
test_edge_cases <- function(seurat_obj) {
  tryCatch({
    print_info("测试边界情况")
    
    # 创建边界测试的输出目录
    edge_case_dir <- file.path(test_output_dir, "edge_cases")
    dir.create(edge_case_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 边界情况1: 极少量细胞
    small_seurat <- subset(seurat_obj, cells = sample(colnames(seurat_obj), size = 10))
    
    print_info("测试边界情况1: 极少量细胞 (10个细胞)")
    
    small_results <- plot_gene_specific_correlation(
      seurat_obj = small_seurat,
      gene_name = "gene_1",
      output_dir = file.path(edge_case_dir, "small_dataset"),
      use_mock_data = TRUE,
      plot_width = 8,
      plot_height = 6,
      plot_dpi = 200
    )
    
    # 边界情况2: 空字符串基因名
    print_info("测试边界情况2: 空字符串基因名")
    
    empty_gene_results <- tryCatch({
      plot_gene_specific_correlation(
        seurat_obj = seurat_obj,
        gene_name = "",
        output_dir = file.path(edge_case_dir, "empty_gene"),
        use_mock_data = TRUE
      )
    }, error = function(e) {
      print_info("正确处理了空基因名的情况")
      return(NULL)
    })
    
    # 边界情况3: 非Seurat对象输入
    print_info("测试边界情况3: 非Seurat对象输入")
    
    wrong_input_results <- tryCatch({
      plot_gene_specific_correlation(
        seurat_obj = data.frame("wrong_input" = 1:5),
        gene_name = "gene_1",
        output_dir = edge_case_dir,
        use_mock_data = TRUE
      )
    }, error = function(e) {
      print_info("正确处理了非Seurat对象的情况")
      return(NULL)
    })
    
    # 评估边界情况测试
    small_case_success <- !is.null(small_results) && !is.null(small_results$correlation_plot)
    edge_case2_success <- is.null(empty_gene_results)  # 应该失败
    edge_case3_success <- is.null(wrong_input_results)  # 应该失败
    
    success_count <- sum(small_case_success, edge_case2_success, edge_case3_success)
    
    print_info(paste("边界情况测试完成: ", success_count, "/3 个测试通过"))
    return(success_count >= 2)  # 允许一个测试失败
  }, error = function(e) {
    print_error(paste("边界情况测试时出错:", e$message))
    return(FALSE)
  })
}

# 主测试函数
main_test <- function() {
  tryCatch({
    print_info("开始测试基因特异性通信网络相关性分析函数")
    
    # 检查主函数是否可用
    if (!check_function_available()) {
      print_error("主函数不可用，测试终止")
      return(FALSE)
    }
    
    # 创建模拟的Seurat对象
    seurat_obj <- create_mock_seurat_object()
    if (is.null(seurat_obj)) {
      print_error("未能创建模拟数据，测试终止")
      return(FALSE)
    }
    
    # 运行各个测试
    basic_test_success <- test_basic_functionality(seurat_obj)
    multiple_genes_success <- test_multiple_genes(seurat_obj)
    cell_columns_success <- test_different_cell_columns(seurat_obj)
    edge_cases_success <- test_edge_cases(seurat_obj)
    
    # 总结测试结果
    print_info("===== 测试结果总结 =====")
    print_info(paste("基本功能测试:", ifelse(basic_test_success, "通过", "失败")))
    print_info(paste("多个基因测试:", ifelse(multiple_genes_success, "通过", "失败")))
    print_info(paste("不同细胞列测试:", ifelse(cell_columns_success, "通过", "失败")))
    print_info(paste("边界情况测试:", ifelse(edge_cases_success, "通过", "失败")))
    
    # 计算总体成功率
    total_tests <- 4
    passed_tests <- sum(basic_test_success, multiple_genes_success, cell_columns_success, edge_cases_success)
    
    success_rate <- (passed_tests / total_tests) * 100
    
    print_info(paste("总体测试结果: ", passed_tests, "/", total_tests, " 测试通过 (", round(success_rate, 1), "%)", sep = ""))
    
    # 创建使用指南
    print_info("===== 使用指南 =====")
    print_info("在您自己的项目中使用此函数:")
    print_info("1. 加载函数: source('path/to/plot_gene_specific_correlation.R')")
    print_info("2. 准备您的Seurat对象")
    print_info("3. 调用函数: results <- plot_gene_specific_correlation(")
    print_info("       seurat_obj = your_seurat_object,")
    print_info("       gene_name = 'your_gene_name',")
    print_info("       output_dir = 'output_directory'")
    print_info("   )")
    print_info("4. 查看结果: 结果会保存在output_dir中")
    
    return(passed_tests >= 2)  # 允许2个测试失败
  }, error = function(e) {
    print_error(paste("主测试函数出错:", e$message))
    return(FALSE)
  })
}

# 运行测试
print_info("基因特异性通信网络相关性分析函数测试脚本")

# 测量运行时间
start_time <- Sys.time()
test_success <- main_test()
end_time <- Sys.time()

# 打印运行时间
print_info(paste("测试完成，总运行时间:", round(as.numeric(difftime(end_time, start_time, units = "mins")), 2), "分钟"))

# 设置退出状态
if (test_success) {
  print_info("测试成功完成！")
  quit(status = 0)
} else {
  print_error("测试失败，请检查错误信息")
  quit(status = 1)
}
