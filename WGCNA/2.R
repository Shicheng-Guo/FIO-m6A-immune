## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


clust = cutreeStatic(sampleTree, cutHeight = 65, minSize = 1)
table(clust)
#clust
#0  1 
#2 98 
keepSamples = (clust==1)
dataExpr = dataExpr[keepSamples, ]
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
