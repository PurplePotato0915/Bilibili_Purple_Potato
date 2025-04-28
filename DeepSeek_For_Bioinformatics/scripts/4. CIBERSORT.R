
# 数据输入
dat_rna = read.csv("transcriptome/0.data/TCGA_mrna.csv", row.names = 1)
dat_clinical = read.csv("transcriptome/0.data/TCGA_clinical.csv", row.names = 1)
dat_clinical = dat_clinical[1:2]

table(rownames(dat_clinical) == names(dat_rna))

str(dat_clinical)
str(dat_rna)

"transcriptome/CIBERSORT.R"
"transcriptome/0.data/LM22.txt"

# 加载所需包
library(ggplot2)      # 绘图
library(ggpubr)       # 统计可视化
library(tidyr)        # 数据整理
library(dplyr)        # 数据处理
library(rstatix)      # 统计检验（包含 add_significance）

# 1. 运行CIBERSORT -----------------------------------------------------------
source("transcriptome/CIBERSORT.R")
lm22 <- normalizePath("transcriptome/0.data/LM22.txt")

# expr_matrix <- as.matrix(t(dat_rna))

# 需要将表达矩阵写入临时文件（CIBERSORT的特殊要求）
expr_file <- "transcriptome/1.tmp/temp_expression_matrix.txt"
write.table(dat_rna, file = expr_file, sep = "\t", quote = FALSE)

# 运行CIBERSORT的正确方式（注意参数顺序）
cibersort_result <- CIBERSORT(
  sig_matrix = lm22,          # 签名矩阵文件路径
  mixture_file = expr_file,   # 表达矩阵文件路径（必须写入文件）
  perm = 1000,
  QN = TRUE
)

# 2. 合并临床数据 -----------------------------------------------------------
clinical_df <- dat_clinical

# 确保样本顺序一致
all(rownames(cibersort_result) == rownames(clinical_df)) # 应该返回TRUE

# 合并结果
merged_data <- cbind(clinical_df, cibersort_result)

# 3. 差异分析 --------------------------------------------------------------
# 转换为长格式
long_data <- merged_data %>%
  select(PATIENT_ID, Group, 1:(ncol(merged_data)-3)) %>% # 排除最后三列（P值等）
  pivot_longer(cols = -c(PATIENT_ID, Group), 
               names_to = "CellType", 
               values_to = "Composition")

# 进行统计检验
stat_test <- long_data %>%
  group_by(CellType) %>%
  t_test(Composition ~ Group, ref.group = "Negative") %>% # 可以调整参考组
  adjust_pvalue(method = "fdr") %>%
  add_significance()

# 筛选显著差异的细胞类型（例如FDR < 0.05）
sig_cells <- stat_test %>%
  filter(p.adj < 0.05) %>%
  pull(CellType)

# 4. 可视化 -----------------------------------------------------------------
stat_labels <- stat_test_sig %>%
  mutate(
    x = 1.5,  # 横坐标居中显示
    label = ifelse(p.adj < 0.001, "***",
                   ifelse(p.adj < 0.01, "**",
                          ifelse(p.adj < 0.05, "*", "")))
  ) %>%
  left_join(
    long_data %>%
      group_by(CellType) %>%
      summarise(y_pos = max(Composition) * 1.1),
    by = "CellType"
  )

p <- ggplot(
  data = long_data %>% 
    filter(CellType %in% sig_cells) %>%
    mutate(Group = factor(Group, levels = c("Negative", "Positive"))),
  aes(x = Group, y = Composition, fill = Group)
) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, shape = 21, color = "black") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
  labs(x = "Group", y = "Cell Composition") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90")
  )

p + geom_text(
  data = stat_labels,
  aes(x = x, y = y_pos, label = label),
  size = 6,
  color = "black",
  inherit.aes = FALSE
) 

# 保存图形
ggsave("final_immune_plot.pdf", width = 12, height = 8)





























