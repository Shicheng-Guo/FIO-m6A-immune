
geneTree=hclust(as.dist(dissTOM),method = "average")

sizeGrWindow(12,9)

plot(geneTree,xlab = "",sub = "",main = "Gene clustering on TOM-based dissimilarity",labels = FALSE,hang = 0.04)

minModuleSize=10

dynamicMods=cutreeDynamic(dendro = geneTree,distM = dissTOM,
            deepSplit = 2,pamRespectsDendro = FALSE,
            minClusterSize = minModuleSize
        )
plot(dynamicMods)
table(dynamicMods)

colors<-labels2colors(dynamicMods)

plotDendroAndColors(geneTree, colors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
