

nSelect=400
set.seed(123)

select=sample(nGenes,size = nSelect)

selectTOM<-dissTOM[select,select]

selectTree=hclust(as.dist(selectTOM),method = "average")

selectColors=colors[select]

sizeGrWindow(9,9)

plotDiss=selectTOM^9

diag(plotDiss)=NA

TOMplot(plotDiss,
        selectTree,
        selectColors,
        main="Network heatmap plot(selected genes)")
