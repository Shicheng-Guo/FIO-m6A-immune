
load("m6a.rda")
tcga_expr <- read.table("CHI3L1.txt", row.names = 1, header = T, check.names = F)
tcga_m6a <- tcga_m6a[colnames(tcga_expr),]
index <- rownames(tcga_expr)
y<-as.numeric(tcga_expr)
head(y)
colnames <- colnames(tcga_m6a)
data <- data.frame(colnames)


for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(tcga_m6a[,i]),y,type="spearman")
  data[i,2] <- test$estimate                                            
  data[i,3] <- test$p.value
}


names(data) <- c("symbol","correlation","pvalue")
head(data)


data %>% 
  #filter(pvalue <0.05) %>% 
  ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
  geom_segment(aes(xend=0,yend=symbol)) +
  geom_point(aes(col=pvalue,size=abs(correlation))) +
  scale_colour_gradientn(colours=c( "#7fc97f","#984ea3")) +
  #scale_color_viridis_c(begin = 0.5, end = 1) +
  scale_size_continuous(range =c(2,8))  +
  theme_minimal() +
  ylab(NULL)+
  xlab("correlation CHI3L1")
