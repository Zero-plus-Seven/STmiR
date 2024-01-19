# 设置工作路径
setwd("D:/jupyter projects/workflow")
library(tidyverse)

#读取CCLE mRNA数据
ccle_mrna <- data.table::fread("../CCLE/CCLE_RNAseq_genes_rpkm_20180929.gct")
#读取CCLE miRNA数据
ccle_mirna <- data.table::fread("../CCLE/CCLE_miRNA_20181103.gct")

# 获取所有列名_后面的字段，并去重
ccle_cancers <- unique(sub("^[^_]+_(.*)", "\\1", colnames(ccle_mrna)))
#删除Description字段
ccle_cancers <- ccle_cancers[-1]

for (cancer in ccle_cancers){
  extracted_column1 <- ccle_mrna[,grepl(paste0(cancer,"$"),colnames(ccle_mrna))]
  cancer_mrna <- ccle_mrna[,extracted_column1,with = FALSE]
  # 删除列名中从 "_" 开始的所有字符,保留样本名称
  colnames(cancer_mrna) <- sub("_.*$", "", colnames(cancer_mrna))
  
  extracted_column2 <- ccle_mirna[,grepl(paste0(cancer,"$"),colnames(ccle_mirna))]
  cancer_mirna <- ccle_mirna[,extracted_column2,with = FALSE]
  # 删除列名中从 "_" 开始的所有字符,保留样本名称
  colnames(cancer_mirna) <- sub("_.*$", "", colnames(cancer_mirna))
  
  #取mrna和mirna的共同样本
  common_columns <- intersect(colnames(cancer_mrna),colnames(cancer_mirna))
  if (length(common_columns) == 0) {
    # 如果长度为0，跳过当前循环，进入下一次循环
    next
  }
  cancer_mrna = cancer_mrna[,common_columns,with = F]
  cancer_mirna = cancer_mirna[,common_columns,with = F]
  #对mrna数据进行对数转换
  cancer_mrna = log1p(cancer_mrna)
  cancer_mrna <- cbind(ccle_mrna$"Description",cancer_mrna)

  #对mirna数据进行对数转换
  cancer_mirna = log1p(cancer_mirna)
  cancer_mirna <- cbind(ccle_mirna$"Description",cancer_mirna)
  
  # 将变量名组合到文件名中
  filename <- paste0("./ccle_cancers/mrna/ccle_mrna_", paste(cancer, collapse="_"), ".csv")
  write.csv(cancer_mrna,filename,row.names = FALSE)
  
  filename <- paste0("./ccle_cancers/mirna/ccle_mirna_", paste(cancer, collapse="_"), ".csv")
  write.csv(cancer_mirna,filename,row.names = FALSE)
}
