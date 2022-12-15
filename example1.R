library(glmnet)

# introdution

# glmnet
startTime=Sys.time()
X <- leukemia$x
y <- leukemia$y
colnames(X) <- paste0("X", 1:3571)
set.seed(555)
fit.glmnet <- glmnet(X, y, family = "binomial", alpha = 0.95,lambda = 10^(seq(0, -2, -0.05)))

# figure1 right
plot(fit.glmnet, "lambda", label=TRUE,xaxt='n')
axis(1,log(fit.glmnet[["lambda"]]),log(fit.glmnet[["lambda"]],10))
grid(lwd = 2)
abline(v=log(0.1584893),col='red',lwd=2,lty=2)
# figure1 bottom left
plot(log(rev(fit.glmnet$lambda),10),rev(fit.glmnet$df),xlab = 'log lambda',ylab = 'Number of nonzero coefficients')
grid(lwd = 2)
abline(v=log(0.1584893,10),col='red',lwd=2,lty=2)
# the number of nonzeros is 6 (with lambda=0.2511886)
# the number of nonzeros is 9 (with lambda=0.1584893)

set.seed(555)
out.cv.glmnet <- cv.glmnet(X, y, family = "binomial", type.measure = "class",lambda = 10^(seq(-2, 0, 0.05)), alpha = 0.95, nfolds = 20)
print(out.cv.glmnet$lambda.1se)
# 0.2511886 (seed=517)
# 0.1584893 (seed=555)

# figure1 top left
plot(log(rev(out.cv.glmnet$lambda),10),rev(out.cv.glmnet$cvm),xlab = 'log lambda',ylab = 'Classification Error: 20-fold cv',ylim = c(0,0.4))
grid(lwd = 2)
abline(v=log(0.1584893,10),col='red',lwd=2,lty=2)
lines(log(rev(out.cv.glmnet$lambda),10),rev(out.cv.glmnet$cvup),col='blue')
lines(log(rev(out.cv.glmnet$lambda),10),rev(out.cv.glmnet$cvlo),col='blue')
lines(log(rev(out.cv.glmnet$lambda),10),rev(error),col='orange',lwd=2)

y.glmnet <- c(predict(fit.glmnet, X, s = out.cv.glmnet$lambda.1se,type = "class"))
print(table(true = factor(y), pred = factor(y.glmnet)))
#     pred
# true  0  1
#   0  47  0
#   1   2 23

error<-c()
for (i in 1:41){
  a<-predict(fit.glmnet, X, s = out.cv.glmnet$lambda[i],type = "class")
  error[i]<-length(which(a!=factor(y)))/72
}
# error
# [1] 0.34722222 0.34722222 0.34722222 0.34722222 0.34722222 0.34722222 0.34722222 0.34722222 0.34722222
# [10] 0.34722222 0.13888889 0.09722222 0.08333333 0.05555556 0.05555556 0.02777778 0.02777778 0.02777778
# [19] 0.02777778 0.02777778 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
# [28] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
# [37] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
endTime=Sys.time()
print(endTime-startTime)
# Time difference of 0.5765512 secs

# varbvs
startTime=Sys.time()
library(varbvs)
set.seed(555)
fit.varbvs <- varbvs(X = X, y = y, Z = NULL, family = "binomial")
endTime=Sys.time()
print(endTime-startTime)
# Time difference of 5.977688 secs
print(summary(fit.varbvs))
# Summary of fitted Bayesian variable selection model:
#   family:     binomial   num. hyperparameter settings: 20
# samples:    72         iid variable selection prior: yes
# variables:  3571       fit prior var. of coefs (sa): yes
# covariates: 1          fit approx. factors (eta):    yes
# maximum log-likelihood lower bound: -28.0321
# Hyperparameters: 
#   estimate Pr>0.95             candidate values
# sa          1.75 [1.62,1.82]         NA--NA
# logodds*   -3.27 [-3.55,-2.88]       (-3.55)--(-1.00)
# *See help(varbvs) for details on how to convert between the
# prior log-odds and the prior inclusion probability.
# Selected variables by probability cutoff:
#   >0.10 >0.25 >0.50 >0.75 >0.90 >0.95 
# 1     1     1     1     1     1 
# Top 5 variables by inclusion probability:
#   index variable    prob PVE    coef*  Pr(coef.>0.95)
# 1  3441    X3441 1.00000  NA -3.20371 [-3.858,-2.549]
# 2  1608    X1608 0.00278  NA -0.00191 [-1.238,-0.147]
# 3  2141    X2141 0.00228  NA -0.00154 [-1.234,-0.123]
# 4  3038    X3038 0.00219  NA  0.00150 [+0.118,+1.259]
# 5  2230    X2230 0.00215  NA  0.00141 [+0.116,+1.208]
# *See help(varbvs) about interpreting coefficients in logistic regression.
m <- length(fit.varbvs$logodds) # default number of setting pi is 20.
err <- rep(0, m)
for (i in 1:m) {
  r <- subset(fit.varbvs, logodds == fit.varbvs$logodds[i])
  ypred <- predict(r, X,type = "class")
  err[i] <- mean(y != ypred)
}

logpi<-log(1/(1+10^(-fit.varbvs$logodds)),10)
# error
plot(logpi,err,xlab = 'log pi',ylab = 'Classification Error',ylim = c(0,0.09),type = 'b',lwd=2)
grid(lwd = 2)
# > err
# [1] 0.08333333 0.08333333 0.08333333 0.08333333 0.06944444 0.05555556 0.05555556 0.05555556 0.05555556
# [10] 0.04166667 0.02777778 0.01388889 0.01388889 0.01388889 0.00000000 0.00000000 0.00000000 0.00000000
# [19] 0.00000000 0.00000000

# expected number of included variables
num<-c()
for (i in 1:20){
  num[i]<-sum(fit.varbvs$alpha[,i])
}
# num
# [1]   1.333460   1.454401   1.619378   1.846757   2.152983   2.576296   3.166647   3.991331   5.189581
# [10]   6.944153   9.631213  13.915290  21.008474  32.548426  52.644515  76.409362 109.853165 155.338582
# [19] 215.938271 295.484892
plot(logpi,num,xlab = 'log pi',ylab = 'expected number of included variables',type = 'b',lwd=2)
grid(lwd = 2)
abline(h=10,col='red',lwd=2,lty=2)
# posterior probability of pi
plot(logpi,fit.varbvs$w,xlab = 'log pi',ylab = 'posterior probability of pi',ylim = c(0,0.25),type = 'b',lwd=2)
grid(lwd = 2)

y.varbvs <- predict(fit.varbvs, X,type = 'class')
print(table(true = factor(y), pred = factor(y.varbvs)))
#      pred
# true  0  1
# 0     44 3
# 1     3 22

