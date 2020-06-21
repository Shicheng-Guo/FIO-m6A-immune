geneModuleMembership = as.data.frame(cor(dataExpr, mergedMEs, use = "p"))
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量

modNames = substring(names(mergedMEs), 3)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");



## 只有连续型性状才能只有计算
## 这里把是否属于 outcome 表型这个变量用0,1进行数值化。
outcome = as.data.frame(dataTraits[,5])
names(outcome) = "outcome"
geneTraitSignificance = as.data.frame(cor(dataExpr, outcome, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(outcome), sep="")
names(GSPvalue) = paste("p.GS.", names(outcome), sep="")
