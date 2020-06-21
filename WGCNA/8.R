
MElist=moduleEigengenes(dataExpr,colors = colors)

MEs=MElist$eigengenes

MEDiss=1-cor(MEs)

METree=hclust(as.dist(MEDiss),method = "average")
