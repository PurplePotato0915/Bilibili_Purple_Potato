
# 加载包----
library(Seurat)
library(SeuratData)
library(tidyverse)
library(harmony)
library(SingleR)
library(SeuratDisk)

# 读入数据----
data1 = Read10X("0.data/data_sc_1_GSE291735/")
Seurat_data1 = CreateSeuratObject(data1)

data2 = read.table("0.data/data_sc_2/counts.txt")
Seurat_data2 = CreateSeuratObject(data2, min.cells = 3)

data3 = Read10X_h5('0.data/data_sc_3/sc_tumor.h5')
Seurat_data3 = CreateSeuratObject(data3)

Convert('0.data/data_sc_4/data.h5ad',dest="h5seurat", assay = "RNA",overwrite = F)
Seurat_data4 = LoadH5Seurat('0.data/data_sc_4/data.h5seurat', meta.data = T)
metadata4 = reticulate::py$metadata





# 下载演示数据----
#InstallData("ifnb")
data("ifnb")
count_ifnb = ifnb@assays$RNA$counts
meta_ifnb = ifnb@meta.data
ifnb = CreateSeuratObject(count_ifnb)
ifnb[['stim']] = meta_ifnb$stim
rm(count_ifnb, meta_ifnb)

# 1. 数据预处理与质控 --------------------------
# 将数据拆分为stimulated和control两组
# ifnb.list <- SplitObject(ifnb, split.by = "stim")


# 2. 整合
###########################  锚点整合  ###########################
# 对每个组进行独立质控
# for (i in 1:length(ifnb.list)) {
#   # 计算线粒体基因百分比
#   ifnb.list[[i]][["percent.mt"]] <- PercentageFeatureSet(ifnb.list[[i]], pattern = "^MT-")
# 
#   # 过滤低质量细胞
#   ifnb.list[[i]] <- subset(ifnb.list[[i]],
#                            subset = nFeature_RNA > 200 &
#                              nFeature_RNA < 2500 &
#                              percent.mt < 5)
# 
#   # 数据归一化
#   ifnb.list[[i]] <- NormalizeData(ifnb.list[[i]], verbose = FALSE)
# 
#   # 寻找高变基因
#   ifnb.list[[i]] <- FindVariableFeatures(ifnb.list[[i]], selection.method = "vst", nfeatures = 2000)
# }
# features <- SelectIntegrationFeatures(object.list = ifnb.list) 
# # 寻找锚点
# immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list,
#                                          anchor.features = features)
# 
# # 整合数据
# combined <- IntegrateData(anchorset = immune.anchors)
#####################################################################


###########################  harmony整合  ###########################
# 选择整合特征基因
ifnb = NormalizeData(ifnb) %>% FindVariableFeatures(nfeatures=3000) %>% ScaleData() %>% RunPCA(verbose=F) 

ifnb = RunHarmony(ifnb, group.by.vars='stim',assay.use='RNA',max.iter.harmony = 20)
ElbowPlot(ifnb,ndims=50,reduction='harmony')
#####################################################################

# 3. 降维与聚类 --------------------------
ifnb = RunTSNE(ifnb, reduction = 'harmony', dims = 1:20) %>% RunUMAP(reduction = 'harmony',dims = 1:20)
ifnb = FindNeighbors(ifnb, dims = 1:20, reduction = 'harmony') %>% FindClusters(resolution = 0.3)
UMAPPlot(ifnb)

# 4.细胞注释
############################  手动注释  ############################
gene.markers.TBNM=c(
  "CD3D","CD3E","CD3G", #T cell
  "CD79A","CD79B","MS4A1", #B cell
  "NKG7","GZMB","GNLY", #NK cell
  "CST3","LYZ","CD68","FCGR3A" #Myeloid cell
)

set.seed(2024)
subobj=subset(ifnb, downsample = 300)
DoHeatmap(subobj,features =gene.markers.TBNM)+NoLegend()
####################################################################


############################  singleR  ############################
ref1 = celldex::MonacoImmuneData()
pred.scRNA = SingleR(test = ifnb@assays$RNA$scale.data, ref = ref1, labels = ref1$label.fine, 
                   clusters = ifnb@active.ident)
new.cluster.ids = pred.scRNA$pruned.labels
names(new.cluster.ids) = levels(ifnb)
ifnb = RenameIdents(ifnb, new.cluster.ids)
ifnb[["cell_type"]] = Idents(ifnb)
####################################################################

UMAPPlot(ifnb)






# 5. 差异细胞比例分析
set.seed(2024)
ID = c(sample(1:5, 6548, replace = T), 
       sample(6:10, 7451, replace = T))
data_analysis = data.frame(as.character(Idents(ifnb)), ifnb$stim, as.character(ID))
names(data_analysis) = c("cell", "intervention", "ID")

# 计算每个患者-细胞类型的计数
cell_proportions <- data_analysis %>%
  # 计算原始计数
  group_by(ID, intervention, cell) %>%
  summarise(count = n(), .groups = 'drop') %>%
  # 计算总细胞数
  group_by(ID, intervention) %>%
  mutate(total = sum(count)) %>%
  # 计算比例
  ungroup() %>%
  mutate(proportion = count / total) %>%
  # 移除总数列
  select(-total)

write.csv(cell_proportions, "1.tmp/cell_proportions.csv", row.names = F)

# 执行统计检验
test_results <- cell_proportions %>%
  # 按细胞类型分组
  group_by(cell) %>%
  # 执行Wilcoxon秩和检验
  summarise(
    p_value = tryCatch({
      wilcox.test(proportion ~ intervention, exact = FALSE)$p.value
    }, error = function(e) NA_real_)) %>%
  # 移除NA值
  filter(!is.na(p_value)) %>%
  # 多重检验校正
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    bonferroni = p.adjust(p_value, method = "bonferroni")
  ) %>%
  # 按FDR值排序
  arrange(fdr)


# 6. 差异表达分析 --------------------------
# 设置默认分组
ifnb_mono = subset(ifnb, subset = cell_type == "Classical monocytes")

Idents(ifnb_mono) = ifnb_mono$stim

cluster.markers <- FindAllMarkers(ifnb_mono, 
                                  only.pos = TRUE,
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25,
                                  test.use = "wilcox")








