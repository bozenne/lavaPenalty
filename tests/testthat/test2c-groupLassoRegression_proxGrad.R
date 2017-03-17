### test-groupLassoRegression_proxGrad.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar  7 2017 (15:44) 
## Version: 
## last-updated: mar 15 2017 (12:33) 
##           By: Brice Ozenne
##     Update #: 33
#----------------------------------------------------------------------
## 
### Commentary: 
##
## BEFORE RUNNING THE FILE:
## library(butils.base) ; package.source("lava.penalty") ;
## path <- file.path(butils.base::path_gitHub(),"lava.penalty","tests")
## source(file.path(path,"FCT.R"))
##
## CONTENT
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(penalized)
library(lava.penalty)
library(testthat)

context("#### Estimate group lasso regression with proximal gradient #### \n")

# {{{ simulation and LVM
# parametrisation
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( c(rep(0,2),1:3) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
categorical(lvm.modelSim, K = 2, labels = letters[1:2]) <- ~X1
categorical(lvm.modelSim, K = 4, labels = letters[1:4]) <- ~X5

# simulation
df.data <- sim(lvm.modelSim,n)
df.data[,c("Y","X2","X3","X4")] <- as.data.frame(scale(df.data[,c("Y","X2","X3","X4")]))

# estimation
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
# }}}

# {{{ no penalty
test_that("NR vs proxGrad with lasso - lambda=0", {
    plvm.model <- penalize(lvm.model)
    
    # likelihood
    eplvm.0 <- estimate(plvm.model, df.data, lambda1 = 0, lambda2 = 0, lambdaG = 0)
    expect_equal(object=coef(elvm.model),
                 expected=coef(eplvm.0),
                 tolerance=test.tolerance, scale=test.scale)    

    # Least square
    eplvm.0 <- estimate(plvm.model,
                        df.data,
                        lambda1 = 0,
                        constrain.lambda = TRUE)
               
    expect_equal(object=coef(elvm.model),
                 expected=coef(eplvm.0),
                 tolerance=test.tolerance,
                 scale=test.scale)    
})
# }}}



#### data ####
library(gglasso)
data(bardet)
group1 <- rep(1, times = 5)#rep(1:5,each=1)
bardet$x <- bardet$x[,1:5]
df.bardet <- data.frame(scale(data.frame(bardet)))

#### models ####

formula_bardety <- as.formula(paste0("y ~ ", paste(names(df.bardet)[names(df.bardet)!="y"], collapse = "+")))
lvm.model_bardety <- lvm(formula_bardety)
lvm.fit_bardety <- estimate(lvm.model_bardety, data = df.bardet)

plvm.model_L1 <- penalize(lvm.model_bardety)
plvm.model_GL <- penalize(lvm.model_bardety) 
plvm.model_GL$penalty$group.coef[] <- 1

gglasso.fit <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls")
seq_lambda <- gglasso.fit$lambda

### no penalization

test_that("LVM vs pLVM with group lasso - lambda=0", {
  plvm.fit_GL <- estimate(plvm.model_GL, data = df.bardet, lambda1 = 0)
  expect_equal(object=coef(plvm.fit_GL),expected=coef(lvm.fit_bardety),tolerance=0.001,scale=NULL)    
})
  

#### group lasso ####
plvm.fit_GL <- estimate(plvm.model_GL,  data = df.bardet, fixSigma = TRUE,
                        lambda1 = 68*seq_lambda[1] * nrow(df.bardet),
                        control = list(constrain = FALSE, iter.max = 1000))
coef(plvm.fit_GL)

plvm.fit_GL <- estimate(plvm.model_GL,  data = df.bardet, fixSigma = TRUE,
                        lambda1 = 69*seq_lambda[1] * nrow(df.bardet),
                        control = list(constrain = FALSE, iter.max = 1000))
coef(plvm.fit_GL)


#### OLD
test <- FALSE
if(test){

mTEST <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls", lambda = 0)
coef(mTEST)

eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet, fixSigma = TRUE,
                           lambda1 = mTEST$lambda[1] * nrow(df.bardet))




iterLambda <- 100
eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet, fixSigma = TRUE,
                           lambda1 = m1$lambda[iterLambda] * nrow(df.bardet),
                           control = list(constrain = FALSE, iter.max = 1000))
coef(eplvm.model_GL)[grep("X",names(coef(eplvm.model_GL)), fixed = TRUE)] - m1$beta[,iterLambda]



#### agreement between lasso and grouped lasso when dealing with one parameter
plvm.model <- penalize(lvm.model_bardety, value = "bardety~X1")

eplvm.model <- estimate(plvm.model,  data = df.bardet,
                           lambda1 = 5,
                           control = list(constrain = TRUE, iter.max = 1000, trace = TRUE))
plvm.model

plvm.model_GL <- plvm.model
plvm.model_GL$penalty$group.penaltyCoef[] <- 1
eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet,
                           lambda1 = 5,
                           control = list(constrain = TRUE, iter.max = 1000, trace = TRUE))
coef(eplvm.model) - coef(eplvm.model_GL)

m1 <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls")

m1$lambda * length(bardet$y) 

m1$b0

m1$beta
m1$npasses

set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)
lvm.test0 <- estimate(lvm.model,  data = df.data)
plvm.test0 <- estimate(plvm.model,  data = df.data)
coef(lvm.test0)

lvm.test1 <- estimate(lvm.model2,  data = df.bardet)
plvm.test1 <- estimate(penalize(lvm.model2),  data = df.bardet)
coef(lvm.test1)
names(df.bardet)

####
library(grplasso)
data(splice)

## Define a list with the contrasts of the factors
contr <- rep(list("contr.sum"), ncol(splice) - 1)
names(contr) <- names(splice)[-1]

## Fit a logistic model 
fit.splice <- grplasso(y ~ ., data = splice, model = LogReg(),
                       contrasts = contr, center = TRUE, standardize = TRUE)

####
library(grpreg)
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X,y,group,penalty="grLasso")
plot(fit)
}

#----------------------------------------------------------------------
### test-groupLassoRegression_proxGrad.R ends here
