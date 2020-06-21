
adjacency<-adjacency(dataExpr,power=9)

TOM<-TOMsimilarity(adjacency)
dissTOM=1-TOM
hierTOM=hclust(as.dist(dissTOM),method = "average")
adj1_cor<-abs(WGCNA::cor(dataExpr,use = "p"))^9
k<-as.vector(apply(adj1_cor,2,sum,na.rm=T))

sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main = "Check Scale Free topology\n")
