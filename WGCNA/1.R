library(WGCNA)

library(reshape2)
library(stringr)

# 
options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads()

exprMat <- "新建文本文档.txt"


type = "unsigned"


corType = "pearson"

corFnc = ifelse(corType=="pearson", cor, bicor)

maxPOutliers = ifelse(corType=="pearson",1,0.05)


robustY = ifelse(corType=="pearson",T,F)

##导入数据##
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                       quote="", comment="", check.names=F)

dim(dataExpr)

## [1] 3600  134

head(dataExpr)[,1:8]





dataExpr<-as.data.frame(t(dataExpr))
