plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

plot(METree)

MEDissThres=0.2

abline(h=MEDissThres,col="red")

merge_modules=mergeCloseModules(dataExpr,colors,cutHeight = 0.2,verbose = 3)


mergedColors=merge_modules$colors

mergedMEs=merge_modules$newMEs

plotDendroAndColors(geneTree,cbind(colors,mergedColors),
                    c("Dynamics Tree Cut","Merged Dynamic"),
                    dendroLabels = FALSE,hang=0.03,
                    addGuide = TRUE,guideHang = 0.05)
