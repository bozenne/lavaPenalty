### test-lassoRegression_EPSODE.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar 16 2017 (09:18) 
## Version: 
## last-updated: mar 17 2017 (19:36) 
##           By: Brice Ozenne
##     Update #: 10
#----------------------------------------------------------------------
## 
### Commentary: 
##
## BEFORE RUNNING THE FILE:
## library(butils.base) ; package.source("lava.penalty") ;
## path <- file.path(butils.base::path_gitHub(),"lava.penalty","tests")
## source(file.path(path,"FCT.R"))
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(penalized)
library(testthat)

#### > settings ####
test.tolerance <- 1e-4
test.scale <- NULL
lava.penalty.options(trace = FALSE, type = "lava")
lava.penalty.options(trace = FALSE, type = "proxGrad")

context("#### Estimate lasso regression with proximal gradient #### \n")


# {{{ simulation and LVM
# parametrisation
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( c(rep(0,2),1:3) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)

# simulation
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))

# estimation
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)

# path
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = FALSE)
seq_lambda1 <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seq_lambda1Sigma <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))
n.lambda <- length(seq_lambda1)

plvm.model <- penalize(lvm.model, Vridge = FALSE)
# }}}

seq_lambda1
# 397.9726 361.0697 348.9896 336.7934 314.8675 284.0589

# {{{ Lars forward
test_that("LARS-forward (lasso)", {

    pathLARS <- estimate(plvm.model,  data = df.data,
                         fit = NULL, regularizationPath = TRUE)
    getPath(pathLARS)
    plot(pathLARS)

    res <- calcLambda(PathFor_fixed)
    res
    plot(res)

    pathEPSODE <- estimate(plvm.model,  data = df.data,
                           fit = NULL, regularizationPath = TRUE,
                           constrain.lambda = TRUE,
                           control.EPSODE = list("resolution_lambda1"=c(0.001,1e-3))
                           )
    getPath(pathLARS, lambda = "lambda1")
    getPath(pathEPSODE, lambda = "lambda1")
    pathEPSODE
    args(getPath.regPath)

    estimate(lvm(Y~X2+X3+X4+X5), data = df.data)
    getPath(PathFor_fixed)
    
    getPath(res, path.constrain = TRUE)
    penalty(res, type = NULL)


    coefType(res)
    
    getPath(PathFor_fixed, lambda = c("lambda1","lambda1.abs"),
            path.constrain = TRUE)

    PathFor_fixed <- estimate(plvm.model,  data = df.data,
                              control = list(constrain = FALSE),
                              regularizationPath = TRUE)

gof(estimate(lvm(Y~X1+X2+X3+X4+X5), data = df.data))$BIC
gof(estimate(lvm(Y~X3+X4+X5), data = df.data))$BIC

 plot(res)
plot(PathFor_fixed)
    plot(PathFor_fixed)
    
    getPath(PathFor_fixed, path.constrain = TRUE)
    getPath(PathFor_fixed, coef = "Y")
    PathFor_fixed
    print(PathFor_fixed, path.constrain = TRUE)
    unclass(x)
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

# {{{ Lars backward

# }}}

# {{{ EPSODE


#----------------------------------------------------------------------
### test-lassoRegression_EPSODE.R ends here
