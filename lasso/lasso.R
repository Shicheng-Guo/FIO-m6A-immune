
library("glmnet")
library("survival")
load("lasso.rda")
cvfit = cv.glmnet(t(myexpr), Surv(mysurv$OS,mysurv$status), 
                  #nfold=10,
                  family = "cox") 
plot(cvfit)
fit <- glmnet(t(myexpr),Surv(mysurv$OS,mysurv$status), 
              family = "cox") 
plot(fit, label = TRUE)
lambda.min=cvfit$lambda.min
coef.min = coef(cvfit, s = "lambda.min") 
coef.min

y=mysurv$status
lasso.prob <- predict(cvfit, newx=t(myexpr) , s=cvfit$lambda.min)
re=cbind(y ,lasso.prob)
head(re)
library(ROCR)
library(caret)

#min
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="blue") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))



###test
################################################################
clinical
exp_data=read.csv("cgga.csv")
exp_data1=read.table("exprcgga.txt",header = T,row.names = 1)
exp_data1=exp_data1[c(c("TAGLN2"  ,  "PDPN"   ,   "TIMP1"   ,  "EMP3","CHI3L1")),rownames(clinical)]
exp_data1=log2(exp_data1+1)
y_test=clinical$Censor
lasso.prob_test <- predict(cvfit, newx=t(exp_data1) , s=cvfit$lambda.min)
re_test=cbind(y_test ,lasso.prob_test)
head(re_test)
library(ROCR)
library(caret)

#min
pred_min <- prediction(re_test[,2], re_test[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="blue") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))




tcga=read.csv("out.csv",header = T,row.names = 1)

q1<-quantile(tcga$riskScore, 0.01)        #取得时1%时的变量值
q99<-quantile(tcga$riskScore, 0.99)       #replacement has 1 row, data has 0 说明一个没换
tcga[tcga$riskScore<q1,]$riskScore<-q1
tcga[tcga$riskScore>q99,]$riskScore<-q99
range(tcga$riskScore)
range(lasso.prob_test)

data=tcga
bestvars<-c("TAGLN2","PDPN","TIMP1","EMP3")
# risk score
rs <- data$riskScore
names(rs) <- rownames(data)
rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
# cutoff  median
rs_data$Risk <- ifelse(rs_data$rs>=median(rs_data$rs), "High-risk", "Low-risk")
head(rs_data)
# follow-up
surv_data <- data.frame(x=1:length(rs),
                        t=data[names(sort(rs)),'OS']*12,
                        s=data[names(sort(rs)),'status']) 
surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
head(surv_data)



exp_data <- data[names(sort(rs)),which(colnames(data) %in% bestvars)]
exp_data[1:2,1:4]


plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=0.5)+
  scale_color_manual(labels=c("High-risk","Low-risk"), 
                     #guide_legend(guide = NULL), 
                     name="Risk score", values =c("#DC0000FF", "#00A087FF")) + 
  
  # 画竖向虚线
  geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                   y = 0, 
                   xend = sum(rs_data$Risk=="Low-risk"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+
  # 画横线
  #geom_segment(aes(x=0,y=median(rs_data$rs),
  #                 xend=nrow(rs_data),
  #                 yend=median(rs_data$rs)),linetype="dashed", size = 0.3)+
  
  # 写文字Cutoff:
  #geom_text(aes(x=sum(rs_data$Risk=="Low-risk")/2,
  #              y=median(rs_data$rs)+8,
  #              label=paste0("Cutoff: ",round(median(rs_data$rs),3))),
  #          col ="black",size = 4,alpha=0.8)+
  
theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  labs(y="Risk score",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
        axis.text.x=element_blank())

plot.A






plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Status),size=0.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")),size=0.6,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels=c("Alive","Dead"),
                     values =c("#00A087FF","#DC0000FF"))+
  labs(y="RFS(months)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
        axis.text.x=element_blank())

plot.B




tmp <- t(scale(exp_data))

tmp<-as.data.frame(tmp)
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp,ifelse(rs_data$Risk=="Low-risk",-1.5,1.5))
A<-rownames(tmp1)
A[5]<-""
tmp2<-cbind(tmp1,A)
tmp2$A
tmp.m <- melt(tmp2,id="A")

p2 <-ggplot(tmp.m, aes(variable, A),size=0.5) + 
  geom_tile(aes(fill = value),position = "identity") 

plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#00A087FF", high="#DC0000FF", mid="white") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())
plot.C


library(patchwork)


plot.A+plot.B+plot.C+plot_layout(ncol = 1)
