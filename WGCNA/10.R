dataTraits<-read.table("c.txt",row.names = 1,header = T)
d<-read.table("mac.txt",header = T,row.names = 1)


rownames(d)<-gsub("\\.","-",rownames(d))
dataTraits<-cbind(dataTraits,d)
dataTraits<-dataTraits[,2:6]
write.csv(dataTraits,"cc.csv")
rownames(MEs)
library(WGCNA)


dataTraits<-dataTraits[rownames(MEs),]
traitColors=numbers2colors(dataTraits,signed = TRUE,centered = TRUE)


plotDendroAndColors(sampleTree,traitColors,groupLabels = names(dataTraits),
                    rowTextAlignment = "right-justified",
                    addTextGuide = FALSE,
                    hang=0.03,
                    dendroLabels = NULL,
                    addGuide = FALSE,
                    guideHang = 0.05,
                    main="Sample dendrogram and trait heatmap")

library(ggplot2)
ggsave("1028.pdf",height = 8,width = 12)
dev.off()
