
##################################  limma  ##################################
library(tidyverse)

dat_rna = read.csv("transcriptome/0.data/TCGA_mrna.csv", row.names = 1)
dat_clinical = read.csv("transcriptome/0.data/TCGA_clinical.csv", row.names = 1)
dat_clinical = dat_clinical[1:2]

str(dat_clinical)
str(dat_rna)

# 加载必要的包
library(limma)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(gridExtra)

# 数据预处理 ----------------------------------------------------------------
# 统一样本ID格式（将临床数据中的"-"替换为"."）
dat_clinical$PATIENT_ID <- gsub("-", ".", dat_clinical$PATIENT_ID)

# 取两个数据集的交集样本
common_samples <- intersect(dat_clinical$PATIENT_ID, colnames(dat_rna))
dat_clinical <- dat_clinical[dat_clinical$PATIENT_ID %in% common_samples, ]
dat_rna <- dat_rna[, common_samples]

# 确保样本顺序一致
dat_clinical <- dat_clinical[match(colnames(dat_rna), dat_clinical$PATIENT_ID), ]

# 创建分组因子（设置Negative为参考组）
Group <- factor(dat_clinical$Group, levels = c("Negative", "Positive"))

# limma差异分析 -----------------------------------------------------------
# 创建设计矩阵
design <- model.matrix(~ Group)

# 拟合线性模型
fit <- lmFit(dat_rna, design)
fit <- eBayes(fit)

# 提取结果
all_genes <- topTable(fit, coef = "GroupPositive", number = Inf, adjust.method = "BH")

# 添加显著性标记（以|logFC|>1 & adj.P.Val<0.05为标准）
all_genes$Significant <- ifelse(
  abs(all_genes$logFC) > 1 & all_genes$adj.P.Val < 0.05,
  "Significant", "Not"
)

# 可视化 -----------------------------------------------------------------
# 1. 火山图
volcano_plot <- ggplot(all_genes, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(adj.P.Value)") +
  theme_bw()

print(volcano_plot)

# 提取前50个显著基因
top50_genes <- rownames(all_genes)[order(all_genes$adj.P.Val)][1:50]
expr_matrix <- as.matrix(dat_rna[top50_genes, ])

# 按照 Group 排序样本
group_order <- order(Group)  # 获取排序索引
expr_matrix <- expr_matrix[, group_order]  # 对表达矩阵按 Group 排序
Group_sorted <- Group[group_order]  # 对 Group 排序

# 绘制热图
pheatmap(
  expr_matrix,
  scale = "row",  # 对行（基因）进行标准化
  cluster_cols = FALSE,  # 禁用列（样本）聚类
  cluster_rows = TRUE,  # 启用行（基因）聚类
  show_rownames = FALSE,  # 不显示行名（基因名）
  show_colnames = FALSE,  # 不显示列名（样本名）
  annotation_col = data.frame(Group = Group_sorted, row.names = colnames(expr_matrix)),  # 添加分组注释
  main = "Top 50 Significant Genes Heatmap (Group Ordered)"
)


# 提取前9个显著基因
top9_genes <- rownames(all_genes)[order(all_genes$adj.P.Val)][1:9]

# 创建绘图数据
plot_data <- data.frame(
  Group = rep(Group, each = length(top9_genes)),
  variable = rep(top9_genes, times = ncol(dat_rna)),
  value = as.vector(t(dat_rna[top9_genes, ]))
)

# 检查 plot_data
print(head(plot_data))
print(str(plot_data))

# 绘制箱线图
plots <- lapply(top9_genes, function(gene) {
  # 提取当前基因的数据
  gene_data <- plot_data[plot_data$variable == gene, ]
  
  # 检查数据是否为空
  if (nrow(gene_data) == 0) {
    stop(paste("No data found for gene:", gene))
  }
  
  # 绘制箱线图
  ggplot(gene_data, aes(x = Group, y = value, fill = Group)) +
    geom_boxplot() +
    labs(title = gene, y = "Expression") +
    theme_minimal() +
    theme(legend.position = "none")
})

# 拼接图形
grid.arrange(grobs = plots, nrow = 3, ncol = 3)

################################## DESeq2  ##################################

# 数据导入----
M_count = read.csv("transcriptome/0.data/M_counts.csv")
M_clinical = read.csv("transcriptome/0.data/M_clinical.csv")

str(M_count)

# 加载必要的包
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# 1. 基因ID转换 -----------------------------------------------------------------
# 提取Entrez ID并转换gene symbol
entrez_ids <- M_count$Entrez_Gene_Id
gene_info <- bitr(entrez_ids, 
                  fromType = "ENTREZID",
                  toType = "SYMBOL",
                  OrgDb = "org.Hs.eg.db")

# 合并symbol信息到表达矩阵
M_count$SYMBOL <- gene_info$SYMBOL[match(M_count$Entrez_Gene_Id, gene_info$ENTREZID)]

# 移除未匹配到的基因和原始Entrez列
M_count <- M_count[!is.na(M_count$SYMBOL), ]
M_count$Entrez_Gene_Id <- NULL

# 处理重复的gene symbol（取表达量总和）
M_count_agg <- aggregate(. ~ SYMBOL, data = M_count, sum)
rownames(M_count_agg) <- M_count_agg$SYMBOL
M_count_agg$SYMBOL <- NULL

# 2. 样本ID统一 -----------------------------------------------------------------
# 处理表达矩阵列名（TCGA.RZ.AB0B.01A -> TCGA-RZ-AB0B）
colnames(M_count_agg) <- gsub("\\.", "-", colnames(M_count_agg))
colnames(M_count_agg) <- substr(colnames(M_count_agg), 1, 12)

# 获取共同样本
common_samples <- intersect(colnames(M_count_agg), M_clinical$PATIENT_ID)

# 过滤并排序数据
M_count_filtered <- M_count_agg[, common_samples]
M_clinical_filtered <- M_clinical[M_clinical$PATIENT_ID %in% common_samples, ]
rownames(M_clinical_filtered) <- M_clinical_filtered$PATIENT_ID
M_clinical_filtered <- M_clinical_filtered[common_samples, ]

# #3. 删除全是 0 或低表达的基因----
# # 计算每个基因的总表达量
# gene_counts <- rowSums(M_count_filtered)
# 
# # 过滤低表达基因（保留总表达量 > 10 的基因）
# M_count_filtered <- M_count_filtered[gene_counts > 10, ]


# 3. DESeq2差异分析 -------------------------------------------------------------
# 创建DESeq对象
dds <- DESeqDataSetFromMatrix(countData = as.matrix(M_count_filtered),
                              colData = M_clinical_filtered,
                              design = ~ Group)

# 预处理：过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= ncol(dds)/2
dds <- dds[keep,]

# 设置因子水平（确保对照是第一个水平）
dds$Group <- relevel(dds$Group, ref = "Low")

# 运行差异分析
dds <- DESeq(dds)
res <- results(dds, contrast = c("Group", "High", "Low"))
resOrdered <- res[order(res$padj), ]

# 4. 可视化 ---------------------------------------------------------------------
# 火山图
EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "High vs Low Group",
                subtitle = "Differential Expression")

# 热图（取padj最小的50个基因）
top50 <- head(order(res$padj), 50)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top50, ]

# 数据标准化（Z-score）
mat <- t(scale(t(mat)))

# 注释信息
annotation_col <- data.frame(Group = colData(dds)$Group)
rownames(annotation_col) <- colnames(mat)

# 绘制热图
pheatmap(mat,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 8,
         main = "Top 50 Significant Genes")

# 输出结果
output = as.data.frame(resOrdered@listData)
output = data.frame(resOrdered@rownames, output)
names(output)[1] = "Gene"
write.csv(output, "transcriptome/2.results/DEGs.csv",row.names = F)



