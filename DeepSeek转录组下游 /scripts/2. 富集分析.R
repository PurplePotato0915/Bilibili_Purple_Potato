
###################################### GO ######################################

# 数据输入
DEGs = read.csv("transcriptome/2.results/DEGs.csv")
str(DEGs)


# 安装必要包（若尚未安装）
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ggplot2"))

# 加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

# 步骤1：准备基因列表（假设物种是人类）
# 提取差异基因（示例筛选条件：padj < 0.05）
significant_genes <- subset(DEGs, padj < 0.01)$Gene

# 将Gene Symbol转换为Entrez ID（需要注意基因命名一致性）
gene.df <- bitr(significant_genes, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

# 步骤2：GO富集分析
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "ALL",       # 可选"BP","MF","CC"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# 步骤3：结果处理
# 提取显著富集结果
ego_result <- as.data.frame(ego)
write.csv(ego_result, "transcriptome/2.results/GO_enrichment_result.csv", row.names = FALSE)

# 步骤4：可视化
# 气泡图（展示前20个条目）
dotplot(ego, 
        showCategory = 20, 
        font.size = 10,
        title = "GO Biological Process Enrichment",
        color = "p.adjust",    # 颜色映射参数
        size = "Count") +       # 点大小映射参数
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal()

# 保存图片
ggsave("GO_enrichment_bubble.pdf", width = 10, height = 8)
ggsave("GO_enrichment_bubble.png", width = 10, height = 8, dpi = 300)



###################################### KEGG ######################################



