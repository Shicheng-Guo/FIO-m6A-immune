
library("glmnet")
library("survival")
load("lasso.rda")
cvfit = cv.glmnet(t(myexpr), Surv(mysurv$OS,mysurv$status), 
                  nfold=10,
                  family = "cox") 
plot(cvfit)
fit <- glmnet(t(myexpr),Surv(mysurv$OS,mysurv$status), 
              family = "cox") 
plot(fit, label = TRUE)
lambda.min=cvfit$lambda.min
coef.min = coef(cvfit, s = "lambda.min") 
coef.min
