library(WGCNA)

library(reshape2)
library(stringr)

# 
options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads()

exprMat <- "新建文本文档.txt"


# 官方推荐 "signed" 或 "signed hybrid"
# 为与原文档一致，故未修改 
type = "unsigned"

# 相关性计算
# 官方推荐 biweight mid-correlation & bicor
# corType: pearson or bicor
# 为与原文档一致，故未修改
corType = "pearson"

corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)

##导入数据##
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                       quote="", comment="", check.names=F)

dim(dataExpr)

## [1] 3600  134

head(dataExpr)[,1:8]





dataExpr<-as.data.frame(t(dataExpr))
