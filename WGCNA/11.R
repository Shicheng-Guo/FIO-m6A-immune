


dataTraits<-dataTraits[rownames(dataExpr),]




moduleTraitCor_noFP<-cor(mergedMEs,dataTraits[,1:5],use = "p")

moduleTraitPvalue_noFP=corPvalueStudent(moduleTraitCor_noFP,nSamples)


textMatrix_noFP<-paste(signif(moduleTraitCor_noFP,2),"\n(",signif(moduleTraitPvalue_noFP,1),")",sep = "")

dev.off()
pdf(file="1029.pdf", paper="letter",width = 7,height =10)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor_noFP,
               xLabels = names(dataTraits[,1:5]),
               yLabels = names(mergedMEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_noFP,
               setStdMargins = FALSE,
               cex.text = 0.65,
               zlim=c(-1,1),
               main=paste("Module-Trait Relationships"))
dev.off()
