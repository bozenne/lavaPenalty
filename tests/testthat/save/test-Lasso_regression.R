### test-Lasso_regression.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 10 2017 (14:23) 
## Version: 
## last-updated: mar  7 2017 (16:33) 
##           By: Brice Ozenne
##     Update #: 14
#----------------------------------------------------------------------
## 
### Commentary: 
##%% standard regression
##> no penalty: lvm-plvm   status=ok
##> lasso path: penalized-plvm   status=~ | WARNING: discrepancy if sigma is not fixed
##> proxGrad at the breakpoints: penalized-plvm  status=ok
##> between breakpoints: lvm-plvm   status=ok | NOTE: if sigma is not fixed, a different solution may be found
##%% partial penalization  
##> status=ok
##%% regression with factors
##> status=~ | NOTE: if sigma is not fixed, a different solution may be found
##%% high dimensional case 
##> proxGrad at the breakpoints: penalized-plvm  status=ok | NOTE: if sigma is not fixed, a different solution may be found. Accuracy is lower compared to the low dimensional case
##> lasso path: penalized-plvm   status=~ | WARNING: discrepancy if sigma is not fixed
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

path <- file.path(butils.base::path_gitHub(),"lava.penalty","tests")
path.res <- file.path(path, "Results/Lasso")
# library(butils.base)
# package.source("lava.penalty")

library(penalized)
library(testthat)
source(file.path(path,"FCT.R"))

#### > settings ####
test.tolerance <- 1e-4
test.scale <- NULL
save <- FALSE
lava.options(trace = FALSE)


# {{{ lasso - linear regression 

# {{{ simulation
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( c(rep(0,2),1:3) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)

# original dataset
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))
df.data2 <- df.data

# same with negative coefficients
df.data2$X5 <- -df.data2$X5
df.data2$X3 <- -df.data2$X3
# }}}
# {{{ models (gold standard)
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = FALSE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
n.lambda <- length(seq_lambda)
# }}}

# {{{ EPSODE forward
test_that("EPSODE-forward vs penalize with lasso", {
  
  # data1
  PathFor_fixed <- estimate(plvm.model,  data = df.data, increasing = TRUE, fixSigma = TRUE, fit = NULL,
                            regularizationPath = TRUE)
  
  lambda1path <- getPath(PathFor_fixed, names = "lambda1.abs")[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance, scale=test.scale)    
  
  if(!is.null(save)){
    if(save){
      saveRDS(PathFor_fixed, file.path(path.res, "RegLD-PathFor_fixed.rds"))   
    }else{
      GS <- readRDS(file.path(path.res, "RegLD-PathFor_fixed.rds"))
      expect_equal(getPath(PathFor_fixed), expected=getPath(GS), tolerance=test.tolerance, scale=test.scale)   
    }
  }
  
  # data2
  PathFor_fixed2 <- estimate(plvm.model,  data = df.data2, increasing = TRUE, estimator = "gaussian", fixSigma = TRUE, 
                              regularizationPath = TRUE)
  
  expect_equal(getPath(PathFor_fixed, names = "lambda1.abs"), getPath(PathFor_fixed2, names = "lambda1.abs"), tolerance=test.tolerance, scale=test.scale)
  
  if(!is.null(save)){
    if(save){
      saveRDS(PathFor_fixed2, file.path(path.res, "RegLD-PathFor_fixed2.rds"))   
    }else{
      GS <- readRDS(file.path(path.res, "RegLD-PathFor_fixed2.rds"))
      expect_equal(getPath(PathFor_fixed2), expected=getPath(GS), tolerance=test.tolerance, scale=test.scale)   
    }
  }
  
})
# }}}

# {{{ EPSODE backward
test_that("EPSODE-backward vs penalize with lasso", {
  
  # data1
  PathBack_fixed <- estimate(plvm.model,  data = df.data, increasing = FALSE, fixSigma = TRUE,
                            regularizationPath = TRUE)
  
  lambda1path <- getPath(PathBack_fixed, names = "lambda1.abs")[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance, scale=test.scale)    
  
  if(!is.null(save)){
    if(save){
      saveRDS(PathBack_fixed, file.path(path.res, "RegLD-PathBack_fixed.rds"))   
    }else{
      GS <- readRDS(file.path(path.res, "RegLD-PathBack_fixed.rds"))
      expect_equal(getPath(PathBack_fixed), expected=getPath(GS), tolerance=test.tolerance, scale=test.scale)   
    }
  }
  
  # data2
  PathBack_fixed2 <- estimate(plvm.model,  data = df.data2, increasing = FALSE, fixSigma = TRUE, 
                             regularizationPath = TRUE)
  
  expect_equal(getPath(PathBack_fixed, names = "lambda1.abs"), getPath(PathBack_fixed2, names = "lambda1.abs"))
  
  if(!is.null(save)){
    if(save){
      saveRDS(PathBack_fixed2, file.path(path.res, "RegLD-PathBack_fixed2.rds"))   
    }else{
      GS <- readRDS(file.path(path.res, "RegLD-PathBack_fixed2.rds"))
      expect_equal(getPath(PathBack_fixed2), expected=getPath(GS), tolerance=test.tolerance, scale=test.scale)   
    }
  }
  
})
# }}}

# {{{ EPSODE partially penalized

penalized.PathL1 <- penalized(Y~.,data = df.data, steps = "Park", trace = FALSE, 
                              unpenalized = Y~ X2 + X4,
                              penalized = Y ~ X1 + X3 + X5)
seq_lambda2 <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

plvm.model2 <- penalize(lvm.model, value = Y ~ X1 + X3 + X5)

test_that("EPSODE-forward vs penalize with partial lasso", {
  system.time(
    PathFor_fixed3 <- estimate(plvm.model2,  data = df.data, increasing = TRUE, fit = NULL,
                                regularizationPath = TRUE)
  )
  lambda1path <- getPath(PathFor_fixed3, names = "lambda1.abs")[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda2,'-'), 2, function(x){which.min(abs(x))})
  expect_equal(lambda1path[indexLambda], expected=seq_lambda2, tolerance=test.tolerance, scale=test.scale)  
  
  if(!is.null(save)){
    if(save){
      saveRDS(PathFor_fixed3, file.path(path.res, "RegLD-PathFor_fixed3.rds"))   
    }else{
      GS <- readRDS(file.path(path.res, "RegLD-PathFor_fixed3.rds"))
      expect_equal(getPath(PathFor_fixed3), expected=getPath(GS), tolerance=test.tolerance, scale=test.scale)   
    }
  }
})

test_that("EPSODE-backward vs penalize with partial lasso", {
  system.time(
    PathBack_fixed3 <- estimate(plvm.model2,  data = df.data, increasing = FALSE, fit = NULL,
                               regularizationPath = TRUE)
  )
  lambda1path <- getPath(PathBack_fixed3, names = "lambda1.abs")[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda2,'-'), 2, function(x){which.min(abs(x))})
  expect_equal(lambda1path[indexLambda], expected=seq_lambda2, tolerance=test.tolerance, scale=test.scale)  
  
  if(!is.null(save)){
    if(save){
      saveRDS(PathBack_fixed3, file.path(path.res, "RegLD-PathBack_fixed3.rds"))   
    }else{
      GS <- readRDS(file.path(path.res, "RegLD-PathBack_fixed3.rds"))
      expect_equal(getPath(PathBack_fixed3), expected=getPath(GS), tolerance=test.tolerance, scale=test.scale)   
    }
  }
})

# }}}

#### > regression with factors ####
cat("* regression with factors \n")

set.seed(10)
n <- 300
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:12), collapse = "+")))
mSim <- lvm(formula.lvm)
categorical(mSim,labels=c("A","B","C")) <- "X1"
categorical(mSim,labels=c("A","B","C")) <- "X11"
df.data <- sim(mSim,n)
df.data[,setdiff(names(df.data),c("X1","X11"))] <- scale(df.data[,setdiff(names(df.data),c("X1","X11"))])

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)

# CRASH
# penalized.PathL1.factor <- penalized(Y~.,data = df.data, steps = "Park", trace = TRUE, lambda2 = 0)

## no penalty
cat("  - no penalty \n")
test_that("NR vs proxGrad with lasso (factor) - lambda=0", {
  elvm.model <- estimate(lvm.model,  data = df.data)
  eplvm.model <-  estimate(plvm.model,  data = df.data, lambda1 = 0, control = list(constrain = TRUE))
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model), tolerance=test.tolerance, scale=test.scale)    
  
  eplvm.model <- estimate(plvm.model, df.data, lambda1 = 0, fixSigma = TRUE, control = list(constrain = TRUE))
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model), tolerance=test.tolerance, scale=test.scale)
})

## regularization path 
cat("  - lasso path \n")
test_that("EPSODE with factors - lambda=0", {
  path1F <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fixSigma = TRUE)
  
  path1B <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fixSigma = TRUE, 
                     increasing = FALSE, resolution_lambda1 = c(1,1e-3))
  p1 <- getPath(path1F)
  p2 <- getPath(path1B, row = 1:nrow(p1))
  rownames(p2) <- 1:nrow(p1)
  expect_equal(p1,p2, tolerance=10*test.tolerance, scale=test.scale)
  
  # comparison to proxGrad - an error could occur here as group lasso is used for factors
  pfit_Fixed <- estimate(plvm.model,  data = df.data, 
                         lambda1 = getPath(path1F, names = "lambda1.abs", row = 8), 
                         fixSigma = TRUE, control = list(constrain = TRUE))
  expect_equal(as.double(coef(pfit_Fixed)), as.double(getPath(path1F, lambda = NULL, row = 8)), tolerance=test.tolerance, scale=test.scale)
  
  ## find another solution when free
  # pfit_Free <- estimate(plvm.model,  data = df.data, lambda1 = getPath(path1F, names = "lambda1", row = 8), fixSigma = FALSE, control = list(constrain = TRUE))
  # expect_equal(as.double(coef(pfit_Free)), as.double(getPath(path1F, lambda = NULL, row = 8)), tolerance=test.tolerance, scale=test.scale)
})

#### > high dimensional ####
cat("* high dimensional case \n")

#### simulation ####
set.seed(10)
n <- 20
n.region <- 25
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:n.region), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( rbinom(n.region, size = 1, prob = 0.3)*rnorm(n.region, 3,1) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))

### models ####
lvm.model <- lvm(formula.lvm)
# elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, lambda1 = 1, steps = "Park", trace = FALSE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

#### check fix lambda ####
cat("  - at the breakpoints \n")

for(iter_l in 1:length(seq_lambda)){
  
   # fixed sigma
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE, method.proxGrad = "FISTA_Beck",
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1,
                               control = list(constrain = TRUE))
  
  # print(range(coef(eplvm.fit_tempo2)-coef2.penalized( penalized.PathL1[[iter_l]])))
  test_that("penalized vs pLVM with lasso (high dimensional - sigmaFixed)", {
    expect_equal(object=coef(eplvm.fit_tempo2),
                 expected=coef2.penalized( penalized.PathL1[[iter_l]]),
                 tolerance=1e-2)
  })
  
  # normal model - can find unexpected solution
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                               control = list(constrain = TRUE))
  
  test_that("penalized vs pLVM with lasso (high dimensional - sigmaFree)", {
    penalized.L1bis <- penalized(Y ~  ., data = df.data, lambda1 = eplvm.fit_tempo1$penalty$lambda1.abs, trace = FALSE)
    expect_equal(object=coef(eplvm.fit_tempo1),
                 expected=coef2.penalized( penalized.L1bis ),
                 tolerance=1e-2)
  })
 
}

#### regularization path ####
cat("  - lasso path \n")
test_that("LVM(EPSODE-backward) vs penalize with lasso (high dimensional)", {
  elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, increasing = FALSE, fixSigma = TRUE,
                                 regularizationPath = TRUE, stopLambda = seq_lambda[5])
  lambda1path <- getPath(elvm.PathL1_EPSODE, names = "lambda1.abs")[,1]
  
  #expect_equal(lambda1path[indexLambda], expected=seq_lambda[-length(seq_lambda)], tolerance=10*test.tolerance, scale=test.scale)    
  
})





n <- 10
x <- rnorm(n)

sd(x)/sqrt(n-1)
sqrt(sum((x-mean(x))^2))/(n-1)

# }}}


#----------------------------------------------------------------------
### test-Lasso_regression.R ends here
