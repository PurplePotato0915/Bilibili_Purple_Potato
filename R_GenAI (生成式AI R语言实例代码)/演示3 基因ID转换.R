
library(biomaRt)
library(dplyr)

# 读取CSV文件
data <- read.csv("Ens_to_Gs.csv")

# 连接Ensembl生物信息数据库
mart <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset("hsapiens_gene_ensembl", mart)

# 执行基因名转换
gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = data$Ens,
  mart = ensembl
)

# 合并数据框，添加Gene Symbol列
data <- data %>%
  left_join(gene_symbols, by = c("Ens" = "ensembl_gene_id"))

# 输出结果
write.csv(data, "Ens_to_Gs_converted.csv", row.names = FALSE)




