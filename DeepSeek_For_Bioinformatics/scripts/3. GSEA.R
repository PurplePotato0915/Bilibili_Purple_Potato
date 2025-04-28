
# 数据输入
DEGs = read.csv("transcriptome/2.results/DEGs.csv")
str(DEGs)

Path_gmt = "transcriptome/0.data/h.all.v7.4.symbols.gmt"

# 加载必要的包
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 准备排序基因列表（按log2FoldChange降序排列）
geneList <- DEGs$log2FoldChange
names(geneList) <- DEGs$Gene
geneList <- sort(geneList, decreasing = TRUE)

# 读取基因集文件
gmt <- read.gmt(Path_gmt)

# 进行GSEA分析
set.seed(123)
gsea_res <- GSEA(
  geneList,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE = gmt
)

# 查看简要结果
head(gsea_res@result)

# 保存结果
write.csv(gsea_res@result, "transcriptome/2.results/GSEA_results.csv")

# 可视化 ----------------------------------------------------------

# 1. 点图（展示最显著通路）
dotplot(gsea_res, showCategory=15, split=".sign") + 
  facet_grid(.~.sign) +
  ggtitle("GSEA - Hallmark Gene Sets")

# 2. 山脊图（展示NES分布）
ridgeplot(gsea_res) + 
  ggtitle("Ridge Plot of GSEA Results")

# 3. 单个通路的GSEA图（选择最显著的通路）
# 提取前3个通路ID
top_pathways <- head(gsea_res$ID, 3)

# 绘制GSEA图
gseaplot2(gsea_res, 
          geneSetID = top_pathways,
          pvalue_table = TRUE,
          title = "Top Enriched Pathways")




