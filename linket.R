library(CellChat)
library(linkET)
library(ggplot2)
library(dplyr)
library(Seurat)

# 假设您已经有了一个Seurat对象（seurat_obj）和一个CellChat对象（cellchat）
# 并且知道要分析的基因名称（例如"MY_GENE"）

# 1. 根据目标基因表达水平对细胞进行分组
# 提取目标基因的表达数据
gene_expression <- FetchData(seurat_obj, vars = "rna_DOCK5")

# 计算基因表达的中位数
median_expression <- median(gene_expression$rna_DOCK5)

# 根据中位数将细胞分为高表达和低表达两组
high_expression_cells <- rownames(gene_expression)[gene_expression$rna_DOCK5 > 0]
low_expression_cells <- rownames(gene_expression)[gene_expression$rna_DOCK5 == 0]

# 2. 创建高表达组的CellChat对象
# 提取高表达组的细胞数据
seurat_high <- subset(seurat_obj, cells = high_expression_cells)

# 创建CellChat对象
cellchat_high <- createCellChat(object = seurat_high, group.by = "cell_subtype") # 替换"cell_subtype"为您的细胞亚型列名

# 设置配体-受体数据库
CellChatDB <- CellChatDB.human # 或 CellChatDB.mouse，根据您的物种
cellchat_high@DB <- CellChatDB

# 预处理数据
cellchat_high <- subsetData(cellchat_high)

# 计算通讯概率
cellchat_high <- identifyOverExpressedGenes(cellchat_high)
cellchat_high <- identifyOverExpressedInteractions(cellchat_high)
cellchat_high <- computeCommunProb(cellchat_high)

# 过滤通讯
cellchat_high <- filterCommunication(cellchat_high, min.cells = 3)

# 计算聚合网络
cellchat_high <- computeCommunProbPathway(cellchat_high)
cellchat_high <- aggregateNet(cellchat_high)

# 3. 创建低表达组的CellChat对象
# 提取低表达组的细胞数据
seurat_low <- subset(seurat_obj, cells = low_expression_cells)

# 创建CellChat对象
cellchat_low <- createCellChat(object = seurat_low, group.by = "cell_subtype") # 替换"cell_subtype"为您的细胞亚型列名

# 设置配体-受体数据库
cellchat_low@DB <- CellChatDB

# 预处理数据
cellchat_low <- subsetData(cellchat_low)

# 计算通讯概率
cellchat_low <- identifyOverExpressedGenes(cellchat_low)
cellchat_low <- identifyOverExpressedInteractions(cellchat_low)
cellchat_low <- computeCommunProb(cellchat_low)

# 过滤通讯
cellchat_low <- filterCommunication(cellchat_low, min.cells = 10)

# 计算聚合网络
cellchat_low <- computeCommunProbPathway(cellchat_low)
cellchat_low <- aggregateNet(cellchat_low)

# 4. 提取通讯概率矩阵
net_high <- cellchat_high@net$prob
net_low <- cellchat_low@net$prob

# 计算总体通讯强度（聚合所有LR对）
overall_comm_high <- apply(net_high, c(1, 2), sum)
overall_comm_low <- apply(net_low, c(1, 2), sum)

# 对称化处理函数
symmetrize_matrix <- function(mat) {
  sym_mat <- (mat + t(mat)) / 2
  diag(sym_mat) <- 0 # 将对角线设为0
  return(sym_mat)
}

sym_matrix_high <- symmetrize_matrix(overall_comm_high)
sym_matrix_low <- symmetrize_matrix(overall_comm_low)

cells_high <- rownames(sym_matrix_high)
cells_low <- rownames(sym_matrix_low)

new_matrix <- matrix(
  data = 0, 
  nrow = length(cells_low), 
  ncol = length(cells_low), 
  dimnames = list(cells_low, cells_low)
)

new_matrix[cells_high, cells_high] <- sym_matrix_high[cells_high, cells_high]
sym_matrix_high=new_matrix
mantel <- mantel_test(sym_matrix_high, sym_matrix_low,
) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
mantel$spec="DOCK5"

# 8. 绘制图形
qcorrplot(correlate(sym_matrix_high), type = "upper", diag = FALSE) +
  # 使用 geom_tile() 替代 geom_square()，或者提供完整的美学映射
  geom_tile() +  # 替代方案1：使用 geom_tile()
  # 或者提供完整的美学映射给 geom_square()
  # geom_square(aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5)) +  # 替代方案2
  
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))


