

FilterGenes_spe<-(abs(geneModuleMembership$MMblue)>0.9)&(abs(geneTraitSignificance$GS.outcome)>0.45)


FilterGenes_spe<-(abs(geneModuleMembership$MMblack)>0.7)
table(FilterGenes_spe)

                                                          
length(FilterGenes_spe)
filter<-(abs(geneModuleMembership$MMblack)>0)
blue<-colnames(dataExpr)[FilterGenes_spe]
write.csv(black,"black.csv")
library(WGCNA)
table(FilterGenes_spe)

hubgene<-colnames(dataExpr)[FilterGenes_spe]
hubgene

hub<-as.data.frame(hubgene)

write.csv(hub,"hubgene.csv")

plotNetworkHeatmap(dataExpr,
                   plotGenes = black,
                   networkType = "signed",
                   useTOM = TRUE,
                   power = 9,
                   main = "signed correlations")

black

ggplot2::ggsave("h.pdf",width = 9,height = 9)
dev.off()
