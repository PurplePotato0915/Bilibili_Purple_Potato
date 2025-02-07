
# 1. 构建数据框
patient_data <- data.frame(
  Case_ID = 1:20,
  Entry_Staff = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
  Stage = c("I", "I", "III", "II", "I", "I", "III", "II", "I", "I", "II", "III", "III", "I", "II", "I", "III", "I", "II", "II")
)

# 2. 生成独一无二的研究号
patient_data$PID <- paste0("PP", sprintf("%03d", 1:20))

# 3. 将分期转换为factor
patient_data$Stage <- as.factor(patient_data$Stage)

# 4. 验证转换成功
str(patient_data)

# 5. 统计不同分期的患者人数
stage_counts <- table(patient_data$Stage)
print(stage_counts)

# 6. 进行Fisher精确检验
contingency_table <- table(patient_data$Entry_Staff, patient_data$Stage)
fisher_result <- fisher.test(contingency_table)
print(fisher_result)

# 7. 导出数据框到当前目录
write.csv(patient_data, file = "patient_data.csv", row.names = FALSE)

patient_data$Entry_Staff
patient_data[1:10,]
patient_data[,2]

