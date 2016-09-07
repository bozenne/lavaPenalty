path <- file.path(butils::dir.gitHub(),"lava.penalty","tests")
path.res <- file.path(path, "Results/ElasticNet")

library(penalized)
library(testthat)
library(lava.penalty)
# package.source("lava.penalty", Rcode = TRUE)
source(file.path(path,"FCT.R"))

context("LVM-lasso")

#### > settings ####
test.tolerance <- 1e-4
test.scale <- NULL
save <- FALSE

#### > factor analysis ####

#### Simulations ####
set.seed(10)
n <- 500
formula.lvm1 <- as.formula( paste0("Y1~eta+", paste0("X",1, collapse = "+") ) )
formula.lvm2 <- as.formula( paste0("Y2~eta+", paste0("X",2, collapse = "+") ) )
formula.lvm3 <- as.formula( paste0("Y3~eta+", paste0("X",3, collapse = "+") ) )


lvm.modelSim <- lvm(list(formula.lvm1,
                         formula.lvm2,
                         formula.lvm3))
latent(lvm.modelSim) <- ~eta
df.data <- sim(lvm.modelSim,n)
df.data <- df.data[,names(df.data) != "eta"]

#### models ####
lvm.model <- lvm(list(Y1~eta,Y2~eta,Y3~eta))
latent(lvm.model) <- "eta"
plvm.model <- penalize(lvm.model)
plvm.modelSim <- penalize(lvm.modelSim)

lvm.ext1 <- lvm(list(Y1~eta+X1+X2,Y2~eta+X2+X3,Y3~eta+X3))
latent(lvm.ext1) <- ~eta
lvm.ext2 <- extendModel(lvm.modelSim, type = "all", covariance = TRUE)


#### no penalty ####
test_that("LVM vs pLVM (lambda = 0)", {
  resLVM <- estimate(lvm.modelSim, data = scale(df.data))
  
  resPLVM1 <- estimate(plvm.modelSim, data = df.data, lambda1 = 0)
  expect_equal(coef(resLVM), coef(resPLVM1), tolerance=test.tolerance, scale=test.scale)
  
  resPLVM2 <- estimate(plvm.modelSim, data = df.data, lambda1 = 0, control = list(constrain = TRUE))
  expect_equal(coef(resLVM), coef(resPLVM2), tolerance=test.tolerance, scale=test.scale)  
  
  resPLVM3 <- estimate(plvm.modelSim, data = df.data, lambda1 = 0, fixSigma = TRUE)
  expect_equal(coef(resLVM), coef(resPLVM3), tolerance=test.tolerance, scale=test.scale)
})

#### A given sigma ####
lambda1 <- 67.24619
test_that("LVM vs pLVM (lambda > 0)", {
  resPLVM1 <- estimate(plvm.modelSim, data = df.data, lambda1 = lambda1, lambda2 = 0)
  resPLVM2 <- estimate(plvm.modelSim, data = df.data, lambda1 = lambda1, lambda2 = 0, control = list(constrain = TRUE))
  expect_equal(coef(resPLVM1), coef(resPLVM2), tolerance=test.tolerance, scale=test.scale)  
  
  resPLVM3 <- estimate(plvm.modelSim, data = df.data, lambda1 = resPLVM2$penalty$lambda1.abs, lambda2 = 0, fixSigma = TRUE)
  expect_equal(coef(resPLVM2), coef(resPLVM3), tolerance=test.tolerance, scale=test.scale)
})

#### grid search
res <- calcLambda(model = plvm.modelSim, seq_lambda1 = seq(0,400, length.out = 10), data.fit = df.data, 
                  keep.fit = TRUE, refit.pLVM = TRUE, 
                  trace = TRUE)
res$regPath$pLVM
plot(seq(0,400, length.out = 10),res$criterion)

#### Regularization path ####
plvm.ext1 <- penalize(lvm.ext1, c("Y1~X1","Y1~X2"))
# plvm.Extended <- penalize(lvm.ext, setdiff(coef(lvm.ext), coef(lvm.model)))
# coef(lvm.modelSim) # 


test_that("pLVM EPSODE vs proxGrad", {
  system.time(
    Path_For <- estimate(plvm.ext1,  data = df.data, increasing = TRUE, fit = NULL, estimator = "numDeriveSimple",
                         regularizationPath = TRUE, lambda2 = 0, stopParam = 5, resolution_lambda1 = c(1e-1,1e-2),
                         control = list(trace =5))
  )
  test <- validPath.lvm(Path_For, data = df.data)
  expect_equal( unique(as.double(test$diff0)), 0)
  expect_equal( test$diff.range, c(0,0), tolerance = 1e-1, scale = NULL)
  
  system.time(             
    Path_Back <- estimate(plvm.ext1,  data = df.data, increasing = FALSE, fit = NULL, estimator = "numDeriveSimple",
                          regularizationPath = TRUE, lambda2 = 0, stopParam = 3, resolution_lambda1 = c(1e-1,1e-2),
                          control = list(trace =TRUE))
  )
  test <- validPath.lvm(Path_Back, data = df.data)
  #expect_equal( test$diff.range, c(0,0), tolerance = 1e-3)
  # expect_equal( test$diff.range, c(0,0), tolerance = 1e-3)
})

#### > factor analysis with factors ####
set.seed(10)
n <- 500
formula.lvm1 <- as.formula( paste0("Y1~eta+", paste0("X",1:2, collapse = "+") ) )
formula.lvm2 <- as.formula( paste0("Y2~eta+", paste0("X",2, collapse = "+") ) )
formula.lvm3 <- as.formula( paste0("Y3~eta+", paste0("X",3, collapse = "+") ) )


lvm.modelSim <- lvm(list(formula.lvm1,
                         formula.lvm2,
                         formula.lvm3))
categorical(lvm.modelSim,labels=c("A","B","C")) <- "X1"
categorical(lvm.modelSim,labels=c("A","B","C")) <- "X3"
latent(lvm.modelSim) <- ~eta
df.data <- sim(lvm.modelSim,n)
df.data <- df.data[,names(df.data) != "eta"]
df.data[,setdiff(names(df.data),c("X1","X3"))] <- as.data.frame(scale(df.data[,setdiff(names(df.data),c("X1","X3"))]))

lvm.model <- lvm(list(formula.lvm1,
                      formula.lvm2,
                      formula.lvm3))
latent(lvm.model) <- ~eta
plvm.model <- penalize(lvm.model)


#### no penalty
test_that("LVM vs pLVM (lambda = 0)", {
  
  resLVM <- estimate(lvm.model, data = df.data)
  
  resPLVM1 <- estimate(plvm.model, data = df.data, lambda1 = 0)
  expect_equal(coef(resLVM), coef(resPLVM1), tolerance=test.tolerance, scale=test.scale)
  
  resPLVM2 <- estimate(plvm.model, data = df.data, lambda1 = 0, control = list(constrain = TRUE))
  expect_equal(coef(resLVM), coef(resPLVM2), tolerance=test.tolerance, scale=test.scale)  
  
  resPLVM3 <- estimate(plvm.model, data = df.data, lambda1 = 0, fixSigma = TRUE)
  expect_equal(coef(resLVM), coef(resPLVM3), tolerance=test.tolerance, scale=test.scale)
})

#### penalty
test_that("LVM vs pLVM (lambda > 0)", {
  resPLVM1 <- estimate(plvm.model, data = df.data, lambda1 = 10)
  
  resPLVM2 <- estimate(plvm.model, data = df.data, lambda1 = 10, control = list(constrain = TRUE))
  expect_equal(coef(resPLVM1), coef(resPLVM2), tolerance=test.tolerance, scale=test.scale)  
  
  resPLVM3 <- estimate(plvm.model, data = df.data, lambda1 = resPLVM2$penalty$lambda1.abs, fixSigma = TRUE)
  expect_equal(coef(resPLVM1), coef(resPLVM3), tolerance=test.tolerance, scale=test.scale)
})
