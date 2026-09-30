########
# set up
########
# Question: how does a change in INCOME affect gas consumption?
# Parameter: average derivative of log gas w.r.t. log income (income elasticity)

rm(list=ls())

library("foreign")
library("dplyr")
library("ggplot2")
library("quantreg")

library("MASS")
library("glmnet")
library("grplasso")
library("nnet")
library("randomForest")
library("gglasso")
library("plotrix")
library("gridExtra")
library("SparseM")

# setwd("path/to/your/repo")   # set this if you are not already in the repo folder

#######################
# clean and format data
#######################
source('get_data.R')
source('primitives.R')
source('stage1.R')
source('stage2.R')

spec=1
#1 means Chernozhukov and Semenova style dictionary, with income as the treatment
#2 means the above plus additional interactions

D_LB=0     #each diagonal entry of \hat{D} lower bounded by D_LB
D_add=.2   #each diagonal entry of \hat{D} increased by D_add
max_iter=10 #max iterations in Dantzig selector iteration over estimation and weights

###########
# algorithm
###########
alpha_estimator=1
gamma_estimator=1
bias=0
#alpha_estimator: 0 dantzig, 1 lasso
#gamma_estimator: 0 dantzig, 1 lasso, 2 rf, 3 nn
#(avoid 2: a random forest barely responds to a small income shift)
#bias: 0 debiased (DML), 1 plug-in

for (quintile in 0:5){   # 0 = full sample, 1-5 = income quintiles

  print(paste0('quintile: '))
  print(paste0(quintile))

  data<-get_data(spec,quintile)

  Y=data[[1]]
  X=data[[2]]
  X.up=data[[3]]
  X.down=data[[4]]
  delta=data[[5]]

  test<-get_MNG(Y,X,X.up,X.down,delta)
  M_hat=test[[1]]
  N_hat=test[[2]]

  # dictionary is identity
  n=nrow(X)
  p=ncol(X)

  #p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
  p0=ceiling(p/4)
  if (p>60){
    p0=ceiling(p/40)
  }

  set.seed(1) # for sample splitting

  results<-rrr(Y,X,X.up,X.down,delta,p0,D_LB,D_add,max_iter,alpha_estimator,gamma_estimator,bias)
  printer(results)
  for_tex(results)
}
