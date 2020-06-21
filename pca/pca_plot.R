save(data,cluster,group,file = "pca.rda")
data.pca <- prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
PCA = data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)
color<-ifelse(group=="cluster1","red","blue")
library(ggplot2)
ggplot(data = PCA, aes(PCA1, PCA2,fill=group)) + geom_point(color=color)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())