#' Cytfly - 通讯概率差异分析 (审美定制版)
#' @export
Cytfly <- function(seurat_obj, target_gene, species, 
                   annotation_col = "seurat_clusters", 
                   output_dir = "./Cytfly_output") {
  
  # 1. 环境准备
  suppressPackageStartupMessages({
    library(Seurat)
    library(CellChat)
    library(linkET)
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
  })
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 2. Seurat v5 兼容处理
  if (as.numeric(substr(as.character(seurat_obj@version), 1, 1)) >= 5) {
    message(">>> 检测到 Seurat v5，正在合并图层...")
    seurat_obj <- tryCatch({ SeuratObject::JoinLayers(seurat_obj) }, error = function(e) seurat_obj)
  }
  
  # 3. 提取数据与分组
  expr_matrix <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  if (!(target_gene %in% rownames(expr_matrix))) {
    target_gene_fix <- gsub("_", "-", target_gene)
    if (!(target_gene_fix %in% rownames(expr_matrix))) stop("找不到基因: ", target_gene)
    target_gene <- target_gene_fix
  }
  gene_vec <- as.numeric(expr_matrix[target_gene, ])
  pos_vals <- gene_vec[gene_vec > 0]
  if (length(pos_vals) < 5) stop("表达细胞不足以分组。")
  seurat_obj$gene_group <- ifelse(gene_vec > median(pos_vals), "high", "low")
  
  # 4. 标签清洗 (增强版：支持注释名称，替换空格和符号)
  labels_vec <- seurat_obj[[annotation_col]]
  if (is.data.frame(labels_vec)) labels_vec <- labels_vec[, 1]
  
  # 关键逻辑：处理注释名称中的非法字符，并加前缀防止 CellChat 报错
  clean_vec <- as.character(labels_vec)
  clean_vec <- gsub("[^[:alnum:]]", "_", clean_vec) # 将非字母数字字符替换为下划线
  seurat_obj$cleaned_labels <- paste0("C_", clean_vec)
  
  # 5. CellChat 循环分析
  cc_list <- list()
  for (grp in c("high", "low")) {
    message("\n>>> 正在处理 [", grp, "] 组通讯概率...")
    sub_obj <- subset(seurat_obj, subset = gene_group == grp)
    if (ncol(sub_obj) < 15) next
    
    tryCatch({
      cc <- createCellChat(object = Seurat::GetAssayData(sub_obj, assay = "RNA", layer = "data"), 
                           meta = data.frame(labels = sub_obj$cleaned_labels, row.names = colnames(sub_obj)), 
                           group.by = "labels")
      cc@DB <- if (tolower(species) == "human") CellChatDB.human else CellChatDB.mouse
      cc_list[[grp]] <- cc %>% subsetData() %>% identifyOverExpressedGenes() %>% 
        identifyOverExpressedInteractions() %>% computeCommunProb() %>% 
        filterCommunication(min.cells = 5) %>% computeCommunProbPathway() %>% aggregateNet()
    }, error = function(e) message("!!! ", grp, " 组失败: ", e$message))
  }
  
  # 6. 提取与对齐概率矩阵
  get_w <- function(cc) { if (is.null(cc) || is.null(cc@net$weight)) return(NULL); as.matrix(cc@net$weight) }
  w_h <- get_w(cc_list$high); w_l <- get_w(cc_list$low)
  all_nodes <- sort(unique(seurat_obj$cleaned_labels))
  
  align_m <- function(m, nodes) {
    out <- matrix(0, nrow=length(nodes), ncol=length(nodes), dimnames=list(nodes, nodes))
    if (!is.null(m)) {
      r <- intersect(rownames(m), nodes); c <- intersect(colnames(m), nodes)
      out[r, c] <- m[r, c]
    }
    return(out)
  }
  final_h <- align_m(w_h, all_nodes); final_l <- align_m(w_l, all_nodes)
  
  # 7. Mantel Test (修复 spec 标签为基因名)
  message(">>> 执行 Mantel Test...")
  m_test <- tryCatch({
    linkET::mantel_test(as.data.frame(final_h), as.data.frame(final_l)) %>%
      dplyr::mutate(
        spec = target_gene, # 核心修改：将左侧标签改为基因名
        rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
        pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
        x = all_nodes[1], y = all_nodes[min(2, length(all_nodes))]
      )
  }, error = function(e) { 
    data.frame(spec=target_gene, env="L", r=0, p=1, rd="< 0.2", pd=">= 0.05", x=all_nodes[1], y=all_nodes[1]) 
  })
  
  # 8. 绘图函数 (移除标题，保留你的审美)
  plot_instance <- function(mat, mantel_df, type) {
    if (sum(mat) == 0) return(NULL)
    
    # 方差过滤
    keep <- apply(mat, 2, sd) > 0
    if (sum(keep) < 2) return(NULL)
    sub_mat <- mat[keep, keep]
    
    # 计算相关性
    corr_res <- linkET::correlate(as.data.frame(sub_mat))
    
    # 绘图逻辑
    p <- qcorrplot(corr_res, type = type, diag = FALSE) +
      geom_tile() +
      geom_couple(aes(colour = pd, size = rd), 
                  data = mantel_df, 
                  curvature = nice_curvature()) +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
      scale_size_manual(values = c("< 0.2" = 0.5, "0.2 - 0.4" = 1, ">= 0.4" = 2)) +
      scale_colour_manual(values = RColorBrewer::brewer.pal(3, "Set2")) +
      guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
             colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
             fill = guide_colorbar(title = "Pearson's r", order = 3))
    # 注意：此处已移除 ggtitle(title)
    
    return(p)
  }
  
  p_h <- plot_instance(final_h, m_test, "upper")
  p_l <- plot_instance(final_l, m_test, "lower")
  
  if(!is.null(p_h)) ggsave(file.path(output_dir, "high_strength_final.png"), p_h, width=9, height=8, dpi=300)
  if(!is.null(p_l)) ggsave(file.path(output_dir, "low_strength_final.png"), p_l, width=9, height=8, dpi=300)
  
  message("=== 分析完成，结果存入：", output_dir)
  return(list(mantel = m_test, prob_h = final_h, prob_l = final_l))
}