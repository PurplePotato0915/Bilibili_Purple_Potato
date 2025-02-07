
# 导入数据
data <- read.csv("pCR.csv")

# 检查缺失值
missing_values <- colSums(is.na(data))

# 输出缺失值情况
print(missing_values)



library(dplyr)

# 将outcome转换为因子
data$outcome <- as.factor(data$outcome)

# 计算中位数
median_results <- data %>%
  group_by(outcome) %>%
  summarise(across(starts_with("cell"), median))

# 输出中位数结果
print(median_results)

# 检验每种细胞在不同终点组的差异（Wilcoxon检验）
p_values <- sapply(data[, 3:12], function(cell) {
  wilcox.test(cell ~ outcome, data = data)$p.value
})

# 输出原始p值
print(p_values)

# 对p值进行FDR校正
p_adjusted <- p.adjust(p_values, method = "fdr")

# 输出校正后的p值
print(p_adjusted)

# 将结果整理成数据框
results <- data.frame(
  Cell = colnames(data)[3:12],
  P_Value = p_values,
  P_Adjusted = p_adjusted
)

# 输出结果
print(results)

# 筛选p.adjust < 0.05的结果
significant_results <- results %>%
  filter(P_Adjusted < 0.05)

# 输出显著性结果
print(significant_results)

# 安装并加载必要的包
library(dplyr)
library(tidyr)

# 提取中位数并重命名列
median_values <- median_results %>%
  pivot_longer(cols = starts_with("cell"), 
               names_to = "Cell", 
               values_to = "Median") %>%
  pivot_wider(names_from = outcome, values_from = Median, 
              names_prefix = "Median_")

# 合并显著性结果并添加中位数列
final_results <- significant_results %>%
  left_join(median_values, by = "Cell")

# 输出最终结果
print(final_results)




# 1. 导入数据并检查缺失值
data <- read.csv("pCR.csv")

# 检查缺失值
missing_values <- colSums(is.na(data))
print(missing_values)

# 2. 计算每种细胞在不同 outcome 组的中位数
# 创建一个空的数据框来存储结果
median_results <- data.frame(Cell_Type = character(), pCR_Median = numeric(), nonpCR_Median = numeric(), stringsAsFactors = FALSE)

# 获取细胞种类的列名
cell_columns <- colnames(data)[3:12]

# 使用 for 循环计算每种细胞的中位数
for (cell in cell_columns) {
  pCR_median <- median(data[data$outcome == "pCR", cell], na.rm = TRUE)
  nonpCR_median <- median(data[data$outcome == "nonpCR", cell], na.rm = TRUE)
  # 将结果添加到数据框中
  median_results <- rbind(median_results, data.frame(Cell_Type = cell, pCR_Median = pCR_median, nonpCR_Median = nonpCR_median))
}

# 打印结果
print(median_results)








