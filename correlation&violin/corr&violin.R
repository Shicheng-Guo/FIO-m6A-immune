load("correlation.rda")
library(ggplot2)
library(ggstatsplot)
#violin
l<-list()

for (i in genes) {
  p=ggpubr::ggviolin(violin,x="immune",y=i)
  l[[i]]=p
}

library(cowplot)
plot_grid(plotlist=l[1:12],ncol = 4,labels = paste0("(",letters[1:12],")"))

#correlation

immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]),y,type="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
library(ggplot2)
library(dplyr)

ggplot(data, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8))+
  #调整legen
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))































immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]),y,type="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
library(ggplot2)
library(dplyr)

ggplot(data, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
