path <- file.path(butils::dir.gitHub(),"lava.penalty","tests")

library(testthat)
library(lava.penalty)
# package.source("lava.penalty", Rcode = TRUE)

context("#### Reduce #### \n")


iter.max <- 1
lambda1 <- 10

m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:2)
m <- regression(m,y='y1',x='z'%++%1:5)

#### simul ####
set.seed(10)
d <- sim(m,5e2)
d <- as.data.frame(scale(d))

suppressWarnings(
  start <- coef(estimate(m, data = d, control = list(iter.max = 0)))
)

## models ####
m.red <- reduce(m)
pm <- penalize(m)
pm.red <- reduce(pm)

#### tests ####
test_that("Regression: plvm vs plvm.reduce", {
  
  method <- lava:::gaussian_method.lvm
  suppressWarnings(
    LVM1 <- estimate(m, data = d,  control = list(iter.max = iter.max, start = start, method = method), estimator = "gaussian1")
  )
  suppressWarnings(
    LVM1.red <- estimate(m.red, data = d, control = list(iter.max = iter.max, start = start[coef(m.red)], method = method), estimator = "gaussian1")
  )
  
  # Hred_E <- gaussian1LP_hessian.lvm(x = m.red, data=d, p=start[coef(m.red)]) ; colnames(Hred_E) <- rownames(Hred_E) <- coef(m.red)
  # Hlava <- lava:::gaussian1_hessian.lvm(x = m,  p=start, n = LVM$data$n, mu = LVM$mu, S = LVM$S)
  # expect_equal(unname(Hlava),unname(Hred_E[coef(m),coef(m)]), tolerance = 1e-6)
  
  expect_equal(coef(LVM1),coef(LVM1.red)[names(coef(LVM1))])
  
  method <- lava:::gaussian2_method.lvm
  suppressWarnings( 
    LVM2 <- estimate(m, data = d, control = list(iter.max = iter.max, start = start, method = method), estimator = "gaussian2")
  )
  # gaussian2_method.lvm
  suppressWarnings(
    LVM2.red <- estimate(m.red, data = d, control = list(iter.max = iter.max, start = start[coef(m.red)], method = method), estimator = "gaussian2")
  )
  
  # Hred_num <- gaussian2LP_hessian.lvm(x = m.red, data=d, p=start[coef(m.red)]) ; colnames(Hred_num) <- rownames(Hred_num) <- coef(m.red) ; attr(Hred_num, "grad") <- NULL
  # Hlava <- lava:::gaussian2_hessian.lvm(x = m,  p=start, n = LVM$data$n, mu = LVM$mu, S = LVM$S, data=d)  ; attr(Hlava, "grad") <- NULL
  # expect_equal(unname(Hlava),unname(Hred_num[coef(m),coef(m)]), tolerance = 1e-6)
  
  expect_equal(coef(LVM2),coef(LVM2.red)[names(coef(LVM2))])
})

test_that("Regression: plvm vs plvm.reduce", {
  method <- lava:::gaussian_method.lvm
  suppressWarnings(
    lassoLVM <- estimate(pm, data = d, lambda1 = lambda1, control = list(iter.max = iter.max, start = start, method = method), estimator = "gaussian1")
  )
  suppressWarnings(
    lassoRLVM <- estimate(pm.red, data = d, lambda1 = lambda1, control = list(iter.max = iter.max, start = start[coef(pm.red)], method = method), estimator = "gaussian1")
  )
  
  expect_equal(coef(lassoRLVM),coef(lassoLVM)[names(coef(lassoRLVM))])
  
  method <- lava:::gaussian2_method.lvm
  suppressWarnings(
    lassoLVM <- estimate(pm, data = d, lambda1 = lambda1, control = list(iter.max = iter.max, start = start[coef(pm)], method = method), estimator = "gaussian2")
  )
  suppressWarnings(
    lassoRLVM <- estimate(pm.red, data = d, lambda1 = lambda1, control = list(iter.max = iter.max, start = start[coef(pm.red)], method = method), estimator = "gaussian2")
  )
  
  expect_equal(coef(lassoRLVM),coef(lassoLVM)[names(coef(lassoRLVM))])
})

test_that("Regression (RP): plvm vs plvm.reduce", {
  
  method <- lava:::gaussian_method.lvm
  RP <- estimate(pm, data = d, regularizationPath = TRUE, estimator = "gaussian1")
  RRP <- estimate(pm.red, data = d, regularizationPath = TRUE, estimator = "gaussian1", fit = NULL)
  # gof(RRP) 
  # vars(pm) ; vars(RRP)
  # pm$M ; RRP$model0$M # RRP$model$model0$M
  
  
  expect_equal(getPath(RP),getPath(RRP)[names(coef(RP))])
})


test_that("Regression: nuclear norm regularization", {
  
})

#### timing
test <- FALSE

if(test == TRUE){
  
  method <- lava:::gaussian_method.lvm
  timout <- 12
  
  system.time( 
    LVM <- try(R.utils:::evalWithTimeout(estimate(m, data = d, control = list(start = start, method = method), estimator = "gaussian2"), timeout = timout))
  )
  # gaussian2_method.lvm
  system.time(
    LVM.red <- try(R.utils:::evalWithTimeout(estimate(m.red, data = d, control = list(start = start[coef(m.red)], method = method), estimator = "gaussian2"), timeout = timout))
  )
  
  
  system.time( 
    pLVM <- try(R.utils:::evalWithTimeout(estimate(pm, data = d, lambda1 = lambda1, control = list(start = start, iter.max = 10, method = method), estimator = "gaussian2"), timeout = timout))
  )
  # gaussian2_method.lvm
  system.time( 
    pLVM.red <- try(R.utils:::evalWithTimeout(estimate(pm.red, data = d, lambda1 = lambda1, control = list(start = start[coef(m.red)], iter.max = 10, method = method), estimator = "gaussian2"), timeout = timout))
  )
  
}

