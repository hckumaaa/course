# example 2

load('/Users/kuma/Desktop/CASIA/final/varbvs-master/data/cfw.RData')
library(varbvs)
library(lattice)
pheno_new<-na.omit(pheno[, c("sacwt", "testis")])
geno_new<-geno[row.names(pheno_new),]
startTime=Sys.time()
fit <- varbvs(geno_new, as.matrix(pheno_new[, "sacwt"]), pheno_new[, "testis"],sa = 0.05, logodds = seq(-5, -3, 0.25),family = 'gaussian')
endTime=Sys.time()
print(endTime-startTime)
# Time difference of 44.72423 secs
print(summary(fit))
# Summary of fitted Bayesian variable selection model:
#   family:     gaussian   num. hyperparameter settings: 9
# samples:    993        iid variable selection prior: yes
# variables:  79748      fit prior var. of coefs (sa): no
# covariates: 2          fit residual var. (sigma):    yes
# maximum log-likelihood lower bound: 2428.7093
# Hyperparameters: 
#   estimate Pr>0.95             candidate values
# sigma   0.000389 [0.000379,0.000404] NA--NA
# sa            NA [NA,NA]             0.05--0.05
# logodds*   -3.78 [-4.25,-3.50]       (-5.00)--(-3.00)
# *See help(varbvs) for details on how to convert between the
# prior log-odds and the prior inclusion probability.
# Selected variables by probability cutoff:
# >0.10 >0.25 >0.50 >0.75 >0.90 >0.95 
# 3     3     3     2     2     1 
# Top 5 variables by inclusion probability:
#   index    variable   prob PVE      coef  Pr(coef.>0.95)
# 1 59249   rs6279141 1.0000  NA -0.008059 [-0.010,-0.006]
# 2 24952  rs33217671 0.9351  NA  0.004761 [+0.003,+0.007]
# 3  9203  rs33199318 0.6869  NA  0.004568 [+0.004,+0.009]
# 4 67415  rs52004293 0.0739  NA  0.000255 [+0.002,+0.005]
# 5 44315 rs253722776 0.0707  NA -0.000259 [-0.006,-0.002]

plot.varbvs(fit,vars = c("rs33199318", "rs33217671", "rs6279141"),var.labels = c("rs33199318", "rs33217671", "rs6279141"),groups = map$chr, gap = 1500)
boxplot(fit$fitted.values-pheno_new[, "testis"],main='boxplot of residuals')
#hist(fit$fitted.values-pheno_new[, "testis"],main='histogram of residuals',xlab = '')
par(mfrow=c(3,3))
for (i in 1:9) {
  #title=paste('theta_',i+': QQplot of predictions and true y')
  qqplot(fit$fitted.values[,i],pheno_new[, "testis"],main=paste('theta_',i,': QQplot of predictions and true y'),xlab = '',ylab = '')
  abline(0,1,col='red')
}

startTime=Sys.time()
name_l<-c("rs33199318", "rs33217671", "rs6279141")
fit2 <- varbvs(geno_new[,name_l], as.matrix(pheno_new[, "sacwt"]), pheno_new[, "testis"],sa = 0.05, logodds = seq(-5, -3, 0.25),family = 'gaussian')
endTime=Sys.time()
print(endTime-startTime)
# Time difference of 0.01645613 secs
print(summary(fit2))
# Summary of fitted Bayesian variable selection model:
#   family:     gaussian   num. hyperparameter settings: 9
# samples:    993        iid variable selection prior: yes
# variables:  3          fit prior var. of coefs (sa): no
# covariates: 2          fit residual var. (sigma):    yes
# maximum log-likelihood lower bound: 2434.1016
# Hyperparameters: 
#   estimate Pr>0.95             candidate values
# sigma   0.000407 [0.000407,0.000408] NA--NA
# sa            NA [NA,NA]             0.05--0.05
# logodds*   -3.06 [-3.25,-3.00]       (-5.00)--(-3.00)
# *See help(varbvs) for details on how to convert between the
# prior log-odds and the prior inclusion probability.
# Selected variables by probability cutoff:
#   >0.10 >0.25 >0.50 >0.75 >0.90 >0.95 
# 3     3     3     3     3     2 
# Top 3 variables by inclusion probability:
#   index   variable  prob PVE     coef  Pr(coef.>0.95)
# 1     3  rs6279141 1.000  NA -0.00824 [-0.010,-0.006]
# 2     2 rs33217671 0.991  NA  0.00520 [+0.003,+0.007]
# 3     1 rs33199318 0.939  NA  0.00648 [+0.004,+0.010]
par(mfrow=c(3,3))
for (i in 1:9) {
  #title=paste('theta_',i+': QQplot of predictions and true y')
  qqplot(fit2$fitted.values[,i],pheno_new[, "testis"],main=paste('theta_',i,': QQplot of predictions and true y'),xlab = '',ylab = '')
  abline(0,1,col='red')
}

par(mfrow=c(3,3))
title(main='scatterplot of predictions from full model and selected model')
for (i in 1:9) {
  #title=paste('theta_',i+': QQplot of predictions and true y')
  plot(fit$fitted.values[,i],fit2$fitted.values[,i],main=paste('theta_',i),xlab = 'full model',ylab = 'selected model')
  abline(0,1,col='red')
}


# example 2.b
logodds <- matrix(-3.78,79748,13)
logodds<- matrix(seq(0,3,0.25) - 3.78,79748,13,byrow = TRUE)
startTime=Sys.time()
fit.geno2 <- varbvs(geno_new, as.matrix(pheno_new[, "sacwt"]), pheno_new[, "testis"], family = "gaussian",logodds = logodds,sa = 0.05)
endTime=Sys.time()
print(endTime-startTime)
# Time difference of 6.636595 mins
# Summary of fitted Bayesian variable selection model:
#   family:     gaussian   num. hyperparameter settings: 13
# samples:    993        iid variable selection prior: no
# variables:  79748      fit prior var. of coefs (sa): no
# covariates: 2          fit residual var. (sigma):    yes
# maximum log-likelihood lower bound: 2428.6612
# Hyperparameters: 
#   estimate Pr>0.95             candidate values
# sigma   0.000386 [0.000381,0.00039]  NA--NA
# sa            NA [NA,NA]             0.05--0.05
# *See help(varbvs) for details on how to convert between the
# prior log-odds and the prior inclusion probability.
# Selected variables by probability cutoff:
#   >0.10 >0.25 >0.50 >0.75 >0.90 >0.95 
# 3     3     3     3     2     2 
# Top 5 variables by inclusion probability:
#   index    variable   prob PVE      coef  Pr(coef.>0.95)
# 1 59249   rs6279141 1.0000  NA -0.008035 [-0.010,-0.006]
# 2 24952  rs33217671 0.9586  NA  0.004869 [+0.003,+0.007]
# 3  9203  rs33199318 0.7536  NA  0.005003 [+0.004,+0.009]
# 4 67415  rs52004293 0.0789  NA  0.000269 [+0.002,+0.005]
# 5 44315 rs253722776 0.0754  NA -0.000273 [-0.005,-0.002]
plot.varbvs(fit.geno2,vars = c("rs33199318", "rs33217671", "rs6279141"),var.labels = c("rs33199318", "rs33217671", "rs6279141"),groups = map$chr, gap = 1500)


library(MASS)
A<-matrix(rnorm(10000,0,1),100,100)
D<-matrix(mvrnorm(100,rep(0,100),diag(100)*2),100,100)
d<-diag(D)
xdx1<-diag(t(A)%*%D%*%A)-(t(A)%*%d)^2/sum(d)
xdx2<-diag(t(A)%*%D%*%A)-(t(A)%*%d/sqrt(sum(d)))^2
xdx1==xdx2
max(abs(xdx1-xdx2))
# 2.842171e-14
