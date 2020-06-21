
cor_ADR<-signif(WGCNA::cor(dataTraits,mergedMEs,use = "p",method = "pearson"),5)
p.values<-corPvalueStudent(cor_ADR,nSamples = nrow(dataTraits))


GS1<-as.numeric(WGCNA::cor(dataTraits[,2],dataExpr,use = "p",method = "pearson"))

Genesignificance<-abs(GS1)

ModuleSignificance<-tapply(Genesignificance,mergedColors,mean,na.rm=T)


ModuleSignificance
