
import scanpy as sc
import pandas as pd

annData = sc.read_h5ad("0.data/data_sc_4/data.h5ad")
metadata = annData.obs

# 方案一 导出并使用R重新读取
metadata.to_csv("D:/Response/datasets/GSE210065/scRNA/metadata.csv",index=False)

# 方案二 在R中直接读取py对象
1
