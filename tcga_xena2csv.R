# 设置工作路径
setwd("D:/jupyter projects/workflow")
library(tidyverse)

#读取临床信息
clin <- data.table::fread("../TCGA PANCAN/Cancer bulk RNA-seq/Survival_SupplementalTable_S1_20171025_xena_sp",data.table = F)
cancer_type <-clin %>% 
  select(sample,'cancer type abbreviation') %>% #选取这些列
  rename(Type='cancer type abbreviation')

#读取pancer mrna数据
tcga_mrna <- read_tsv("../TCGA PANCAN/Cancer bulk RNA-seq/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena")

#读取pancer mirna数据
tcga_mirna <- read_tsv("../TCGA PANCAN/Cancer miRNA-seq/pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena")

#用0替换NA值
tcga_mirna[is.na(tcga_mirna)] <- 0
tcga_mrna[is.na(tcga_mrna)] <- 0

# 将小于0的值替换为0
tcga_mrna[tcga_mrna < 0] <- 0
tcga_mirna[tcga_mirna < 0] <- 0
print(any(tcga_mrna<0))
print(any(tcga_mirna<0))

can_types <- unique(cancer_type$Type)

for (cancer in can_types){
  #取出某一特定癌症类型的mirna数据
  indices <- which(cancer_type$Type == cancer)
  sample_values <- cancer_type$sample[indices]
  #mirna和mrna样本对齐
  common_columns1 <- intersect(colnames(tcga_mrna), as.character(sample_values))
  common_columns2 <- intersect(colnames(tcga_mirna), as.character(common_columns1))
  if (length(common_columns2) == 0) {
    # 如果长度为0，跳过当前循环，进入下一次循环
    next
  }
  
  cancer_mrna = tcga_mrna[common_columns2]
  #对数据进行对数转换
  cancer_mrna = log1p(cancer_mrna)
  cancer_mrna <- cbind(tcga_mrna$"sample",cancer_mrna)
  filename <- paste0("./tcga_cancers/mrna/tcga_mrna_", paste(cancer, collapse="_"), ".csv")
  write.csv(cancer_mrna,filename,row.names = FALSE)
  
  cancer_mirna = tcga_mirna[common_columns2]
  
  #对数据进行对数转换
  cancer_mirna = log1p(cancer_mirna)
  cancer_mirna <- cbind(tcga_mirna$"sample",cancer_mirna)
  # 将变量名组合到文件名中
  filename <- paste0("./tcga_cancers/mirna/tcga_mirna_", paste(cancer, collapse="_"), ".csv")
  write.csv(cancer_mirna,filename,row.names = FALSE)
  
  }

