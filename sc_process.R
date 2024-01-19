setwd("D:/jupyter projects/workflow")
library(Seurat)
library(loomR)
library(SeuratDisk)
outdir = "./prostate_data"

sdata <- Read10X_h5(file = "./prostate_data/PRAD_GSE141445_expression.h5")
data_seurat <- CreateSeuratObject(sdata,project = "data_sample") #后面就可以单细胞处理的标准流程啦
# seurat对象转换为loop文件
sdata.loom <- as.loom(x = data_seurat, filename = "./prostate_data/PRAD_GSE141445_expression.loom")
# Always remember to close loom files when done
sdata.loom$close_all()
