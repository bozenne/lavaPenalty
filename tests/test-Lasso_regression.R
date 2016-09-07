path <- file.path(butils::dir.gitHub(),"lava.penalty","tests")
path.res <- file.path(path, "Results/Lasso")

library(penalized)
library(testthat)
library(lava.penalty)
# package.source("lava.penalty", Rcode = TRUE)
source(file.path(path,"FCT.R"))

context("Reg-lasso")

#### > settings ####
test.tolerance <- 1e-4
test.scale <- NULL
save <- FALSE

#### > standard regression ####

#### simulation ####
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

### models ####
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
n.lambda <- length(seq_lambda)

#### no penalty ####
test_that("NR vs proxGrad with lasso - lambda=0", {
  eplvm.model <- estimate(plvm.model, df.data, lambda1 = 0)
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model), tolerance=test.tolerance, scale=test.scale)    
  
  eplvm.model <- estimate(plvm.model, df.data, lambda1 = 0, fixSigma = TRUE)
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model), tolerance=test.tolerance, scale=test.scale)    
})

#### lasso path ####
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


#### check fix lambda - at the breakpoints ####
pb <- utils:::txtProgressBar(max = n.lambda)

for(iter_l in 1:n.lambda){
  
  # normal model
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2, 
                               control = list(trace = -1))
  
  test_that("NR vs proxGrad with lasso", {
    expect_equal(object = coef(eplvm.fit_tempo1),
                 expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance = test.tolerance, scale=1)    
  })
  
  # with constrains
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2, 
                               control = list(constrain = TRUE, trace = -1))
  
  test_that("NR vs proxGrad with lasso (constrain)", {
    expect_equal(object = coef(eplvm.fit_tempo2),
                 expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance = test.tolerance, scale=1)  
  })
  
  # fixed sigma
  eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1, 
                               control = list(trace = -1))
  
  test_that("NR vs proxGrad with lasso (fix sigma)", {
    expect_equal(object = coef(eplvm.fit_tempo3),
                 expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance = test.tolerance, scale=1)  
  })
  
  utils::setTxtProgressBar(pb, iter_l)
}
close(pb)

#### check fix lambda - between breakpoints ####
seq2_lambda <- igraph::running.mean(seq_lambda, binwidth = 2)

pb <- utils:::txtProgressBar(max = length(seq2_lambda))

for(iter_l in 1:length(seq2_lambda)){
  penalized.L1 <- penalized(Y ~  ., data = df.data, lambda1 = seq2_lambda[iter_l], trace = FALSE)
  
  # constrain
  eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                               lambda1 = seq2_lambda[iter_l], 
                               control = list(trace = -1))
  
  test_that("NR vs proxGrad with lasso (fix sigma - between knots)", {
    expect_equal(object = coef(eplvm.fit_tempo3),
                 expected = coef2.penalized(penalized.L1),
                 tolerance = test.tolerance, scale=1)  
  })
  
  # no constrain - can find unexpected solution
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data,
                              lambda1 = seq2_lambda[iter_l]/penalized.L1@nuisance$sigma2, 
                              control = list(trace = -1))
  
  test_that("NR vs proxGrad with lasso (free sigma - between knots)", {
    penalized.L1bis <- penalized(Y ~  ., data = df.data, lambda1 = eplvm.fit_tempo1$penalty$lambda1.abs, trace = FALSE)
    expect_equal(object = coef(eplvm.fit_tempo1),
                 expected = coef2.penalized(penalized.L1bis),
                 tolerance = test.tolerance, scale=1)
    
    
  })
  utils::setTxtProgressBar(pb, iter_l)
  
}
close(pb)

#### partial penalization ####
ppenalized.PathL1 <- penalized(Y~.,data = df.data, steps = "Park", trace = TRUE, 
                              unpenalized = Y~ X2 + X4,
                              penalized = Y ~ X1 + X3 + X5)
seq_lambda2 <- unlist(lapply(ppenalized.PathL1, function(x){x@lambda1}))

plvm.model2 <- penalize(lvm.model, value = Y ~ X1 + X3 + X5)

test_that("EPSODE-forward vs penalize with partial lasso", {
  system.time(
    PathFor_fixed3 <- estimate(plvm.model2,  data = df.data, increasing = TRUE, fit = NULL,
                                regularizationPath = TRUE,
                                control = list(trace =TRUE))
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
                               regularizationPath = TRUE, 
                               control = list(trace =TRUE))
  )
  lambda1path <- getPath(PathFor_fixed3, names = "lambda1.abs")[,1]
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

#### > regression with factors ####
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
test_that("NR vs proxGrad with lasso (factor) - lambda=0", {
  elvm.model <- estimate(lvm.model,  data = df.data)
  eplvm.model <-  estimate(plvm.model,  data = df.data, lambda1 = 0, control = list(constrain = TRUE))
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model), tolerance=test.tolerance, scale=test.scale)    
  
  eplvm.model <- estimate(plvm.model, df.data, lambda1 = 0, fixSigma = TRUE, control = list(constrain = TRUE))
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model), tolerance=test.tolerance, scale=test.scale)
})

## regularization path
test_that("EPSODE with factors - lambda=0", {
  path1F <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fixSigma = TRUE)
  
  path1B <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fixSigma = TRUE, 
                     increasing = FALSE, resolution_lambda1 = c(1,1e-3))
  p1 <- getPath(path1F, order = "lambda1.abs")
  p2 <- getPath(path1B, order = "lambda1.abs", row = 1:nrow(p1))
  rownames(p2) <- 1:nrow(p1)
  expect_equal(p1,p2, tolerance=10*test.tolerance, scale=test.scale)
  
  # comparison to proxGrad
  pfit_Fixed <- estimate(plvm.model,  data = df.data, 
                         lambda1 = getPath(path1F, names = "lambda1.abs", row = 8), 
                         fixSigma = TRUE, control = list(constrain = TRUE))
  expect_equal(as.double(coef(pfit_Fixed)), as.double(getPath(path1F, getLambda = NULL, row = 8)), tolerance=test.tolerance, scale=test.scale)
  
  ## find another solution when free
  # pfit_Free <- estimate(plvm.model,  data = df.data, lambda1 = getPath(path1F, names = "lambda1", row = 8), fixSigma = FALSE, control = list(constrain = TRUE))
  # expect_equal(as.double(coef(pfit_Free)), as.double(getPath(path1F, getLambda = NULL, row = 8)), tolerance=test.tolerance, scale=test.scale)
})

#### > high dimensional ####

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
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, lambda1 = 1, steps = "Park", trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

#### check fix lambda ####
pb <- utils:::txtProgressBar(max = length(seq_lambda))

for(iter_l in 1:length(seq_lambda)){
  
   # fixed sigma
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1,
                               control = list(constrain = TRUE, trace = -1))
  
  # print(range(coef(eplvm.fit_tempo2)-coef2.penalized( penalized.PathL1[[iter_l]])))
  test_that("penalized vs pLVM with lasso (high dimensional - sigmaFixed)", {
    expect_equal(object=coef(eplvm.fit_tempo2),
                 expected=coef2.penalized( penalized.PathL1[[iter_l]]),
                 tolerance=1e-3)
  })
  
  # normal model - can find unexpected solution
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                               control = list(constrain = TRUE, trace = -1))
  
  test_that("penalized vs pLVM with lasso (high dimensional - sigmaFree)", {
    penalized.L1bis <- penalized(Y ~  ., data = df.data, lambda1 = eplvm.fit_tempo1$penalty$lambda1.abs, trace = FALSE)
    expect_equal(object=coef(eplvm.fit_tempo1),
                 expected=coef2.penalized( penalized.L1bis ),
                 tolerance=1e-2)
  })
  utils::setTxtProgressBar(pb, iter_l)
  
}
close(pb)


#### regularization path ####
test_that("LVM(EPSODE-backward) vs penalize with lasso (high dimensional)", {
  elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, increasing = FALSE, fixSigma = TRUE,
                                 regularizationPath = TRUE, stopLambda = seq_lambda[5],
                                 control = list(trace = TRUE))
  lambda1path <- getPath(elvm.PathL1_EPSODE, names = "lambda1.abs", order = "lambda1.abs")[,1]
  
  #expect_equal(lambda1path[indexLambda], expected=seq_lambda[-length(seq_lambda)], tolerance=10*test.tolerance, scale=test.scale)    
  
})



