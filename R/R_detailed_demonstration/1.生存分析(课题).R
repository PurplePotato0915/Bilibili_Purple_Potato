
# 加载包 ----
library(tidyverse)
library(survival)
library(survminer)
library(rms)

# 读取数据 ----
dat_RNA = read.csv(file = "0.data/TCGA_mrna.csv", row.names = 1)
dat_clinical = read.csv("0.data/TCGA_clinical.csv", row.names = 1)

# 检查数据 ----
str(dat_RNA)
str(dat_clinical)

table(dat_clinical$HER2, useNA = "ifany")
table(dat_clinical$ER_STATUS_BY_IHC, useNA = "ifany")
table(dat_clinical$Stage, useNA = "ifany")
table(dat_clinical$SURGICAL_PROCEDURE_FIRST, useNA = "ifany")
table(is.na(dat_clinical$OS_MONTHS))
table(dat_clinical$OS_STATUS, useNA = "ifany")

# 数据预处理 ----
table(rownames(dat_clinical) == colnames(dat_RNA)) #检查样本是否对齐
sum(is.na(dat_RNA)) #查看缺失样本

dat_clinical_analysis = dat_clinical %>%
  select(AGE,Stage,ER_STATUS_BY_IHC,HER2,SURGICAL_PROCEDURE_FIRST,OS_MONTHS, OS_STATUS) %>% 
  rename(Age = AGE,
         ER = ER_STATUS_BY_IHC, 
         Surgery = SURGICAL_PROCEDURE_FIRST) %>% 
  mutate(HER2 = if_else(HER2 == "Positive", 1, 0)) %>%
  mutate(Surgery = if_else(Surgery == "Modified Radical Mastectomy" | Surgery == "Simple Mastectomy", "Mastectomy", 
         Surgery))

# 临床变量单因素cox回归 ----
cox = coxph(Surv(OS_MONTHS, OS_STATUS) ~ Age, data = dat_clinical_analysis); summary(cox)
cox = coxph(Surv(OS_MONTHS, OS_STATUS) ~ Stage, data = dat_clinical_analysis); summary(cox)
cox = coxph(Surv(OS_MONTHS, OS_STATUS) ~ ER, data = dat_clinical_analysis); summary(cox)
cox = coxph(Surv(OS_MONTHS, OS_STATUS) ~ HER2, data = dat_clinical_analysis); summary(cox)
cox = coxph(Surv(OS_MONTHS, OS_STATUS) ~ Surgery, data = dat_clinical_analysis); summary(cox)

km = survfit(Surv(OS_MONTHS, OS_STATUS) ~ Stage, data = dat_clinical_analysis)
ggsurvplot(km, risk.table = T, pval = T)

# 基因表达单因素cox回归 ----
HR = c()
P = c()

time_1 = Sys.time()
for (i in 1:nrow(dat_RNA)) {
  gene = dat_RNA[i,] %>% t() %>% as.data.frame()
  colnames(gene) = "gene"
  gene = gene %>% mutate(gene = if_else(gene > median(gene), "High", "Low")) %>% unlist()
  
  cox = coxph(Surv(dat_clinical_analysis$OS_MONTHS, dat_clinical_analysis$OS_STATUS) ~ gene)
  tmp = summary(cox)
  
  HR[i] = round(tmp$coefficients[2], 2)
  P[i] = round(tmp$coefficients[5], 3)
}
time_2 = Sys.time()
time_diff_for = time_2 - time_1 #Time difference of 4.374822 mins

colnames(result_uni) = c('HR', "P")
result_uni$P.adjusted <- p.adjust(result_uni$P, method = "BH")
#rownames(result_uni) = rownames(dat_RNA)

# 结果导出 ----
write.csv(result_uni, "2.results/result_uni.csv")






# 制作nomo图 ----
result_uni_order = result_uni[order(result_uni$P.adjusted),]

nomodata=data.frame(dat_clinical_analysis$OS_MONTHS,
                    dat_clinical_analysis$OS_STATUS,
                    dat_clinical_analysis$Age,
                    dat_clinical_analysis$Stage,
                    as.numeric(dat_RNA["ATP1A2",]),
                    as.numeric(dat_RNA["FREM1",]),
                    as.numeric(dat_RNA["GPR171",]))
colnames(nomodata)=c('time','state','Age','Stage','ATP1A2',"FREM1",'GPR171')

dd=datadist(nomodata)
options(datadist="dd")

attach(nomodata);psmmodel=psm(Surv(time+0.01,state)~Age+Stage+ATP1A2+FREM1+GPR171,
                              data=nomodata,dist='lognormal');detach(nomodata)
survival=Survival(psmmodel)
survival1=function(x)survival(36,x)
survival2=function(x)survival(60,x)
survival3=function(x)survival(120,x)
nom=nomogram(psmmodel,fun=list(survival1,survival2,survival3),
             fun.at = c(0.05,seq(0.1,0.9,by=0.1),0.95),
             funlabel=c('3-year','5-year','10-year'));plot(nom)

# 校准曲线 !!!注意，本视频仅做技术展示，校准曲线应该用于验证集!!!----
attach(nomodata);cph_model=cph(Surv(time+0.01,state)~Age+Stage+ATP1A2+FREM1+GPR171,
                       data=nomodata,surv=T,x=TRUE, y=TRUE,time.inc=12*5);detach(nomodata)
calibrate_result = calibrate(cph_model,cmethod='KM',method="boot",u=12*5,m=50,B=300)
plot(calibrate_result,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlab="Nomogram-Predicted Probability of 5-Year RFS",
     ylab="Actual 5-Year RFS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))

# ROC !!!注意，本视频仅做技术展示，ROC应该用于验证集!!!----
y = dat_clinical_analysis$OS_STATUS[!is.na(dat_clinical_analysis$Stage)]
y_ = psmmodel$linear.predictors

roc_model=roc(y, y_)
plot(roc_model);auc(roc_model)


# 附加 for循环的改进 - apply族函数 ----
HR_cal = function(x){
  gene = ifelse(x>median(x), "High", "Low")
  cox = coxph(Surv(dat_clinical_analysis$OS_MONTHS, dat_clinical_analysis$OS_STATUS) ~ gene)
  tmp = summary(cox)
  
  a = round(tmp$coefficients[2], 2)
  b = round(tmp$coefficients[5], 3)
  
  return(c(a, b))
}

time_1 = Sys.time()
result_uni = apply(dat_RNA, 1, HR_cal) %>% t()
time_2 = Sys.time()
time_diff_for = time_2 - time_1 #Time difference of 1.163225 mins


# 附加 for循环的进一步改进 - 多线程 ----
library(parallel)
library(future)
library(future.apply)

plan(multisession, workers = detectCores() - 1)

time_1 = Sys.time()
result_uni <- future_apply(dat_RNA, 1, HR_cal) %>% t() %>% as.data.frame()
time_2 = Sys.time()
time_diff_for = time_2 - time_1 #Time difference of 21.07068 secs

plan(sequential)

