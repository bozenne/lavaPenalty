path <- file.path(butils::dir.gitHub(),"lava.penalty","tests")

library(testthat)
library(lava.penalty) # butils:::package.source("lava.penalty", Rcode = TRUE, RorderDescription = FALSE)


context("#### Reduce #### \n")

iter.max <- 1
lambda1 <- 10


#### Regression ####
m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:2)
m <- regression(m,y='y1',x='z'%++%1:5)

## simul
set.seed(10)
d <- sim(m,5e2)
d <- as.data.frame(scale(d))

suppressWarnings(
  start <- coef(estimate(m, data = d, control = list(iter.max = 0)))
)

## models
m.red <- reduce(m)
pm <- penalize(m)
pm.red <- reduce(pm)

## tests moment
test_that("Regression: moment reduce", {
  e <- estimate(m, d, estimator = "gaussian1")
  index <- match(coef(m.red),coef(m))
  
  g1 <- gaussian1LP_gradient.lvm(m.red, p = start[index], data = d, indiv = FALSE)
  g2 <- lava:::gaussian1_gradient.lvm(x = m, data=d, p=start, S = e$S, n = e$data$n, mu = e$mu)

  expect_equal(unname(g1), g2[index])
})

## tests estimation
test_that("Regression: lvm vs lvm.reduce", {
  
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
  RP <- estimate(pm, data = d, regularizationPath = TRUE, estimator = "gaussian1", fit = NULL)
  RRP <- estimate(pm.red, data = d, regularizationPath = TRUE, estimator = "gaussian1", fit = NULL)

  expect_equal(getPath(RP),getPath(RRP)[names(getPath(RP))], tolerance = 1e-6)
})


test_that("Regression: nuclear norm regularization", {
  
})

#### Latent variable model ####
m <- lvm()
m <- regression(m,y=c('y1','y2','y3','y4'),x='eta')
m <- regression(m,y=c('y2','y3'),x='x'%++%1:5)
latent(m) <- ~eta
m <- regression(m,y=c('y4','y2'),x='z'%++%1:2)
covariance(m) <- y2~y1

## simul
set.seed(10)
d <- sim(m,5e2)
d <- as.data.frame(scale(d))

start <- setNames(rep(0, length(coef(m))), coef(m))
start[grep("~", names(start))] <- 1
suppressWarnings(
  startLVM <- coef(estimate(m, data = d, control = list(iter.max = 0)))
)
start[names(startLVM)] <- startLVM

## models
m.red <- reduce(m)
pm <- penalize(m)
pm.red <- reduce(pm)

suppressWarnings(
  start2 <- estimate(m.red, data = d, control = list(iter.max = 0), quick = TRUE)
)

e <- estimate(m, d, estimator = "gaussian1")
index <- match(coef(m.red),coef(m))

## tests moment
test_that("LVM: moment reduce", {
  
 
  g1 <- gaussian1LP_gradient.lvm(m.red, p = start[index], data = d, indiv = FALSE)
  g2 <- lava:::gaussian1_gradient.lvm(x = e$model, data=d, p=startLVM, S = e$S, n = e$data$n, mu = e$mu)
  g2 <- setNames(g2, names(startLVM))
  
  expect_equal(g1[intersect(names(g1),names(g2))], 
               g2[intersect(names(g1),names(g2))])
})

## tests estimation
test_that("LVM: lvm vs lvm.reduce", {
  
  method <- lava:::gaussian_method.lvm
  suppressWarnings(
    LVM1 <- estimate(m, data = d,  control = list(iter.max = iter.max, start = startLVM, method = method), estimator = "gaussian1")
  )
  suppressWarnings(
    LVM1.red <- estimate(m.red, data = d, control = list(iter.max = iter.max, start = start[index], method = method), estimator = "gaussian1")
  )
  
  expect_equal(coef(LVM1)[intersect(names(coef(LVM1)),names(coef(LVM1.red)))],
               coef(LVM1.red)[intersect(names(coef(LVM1)),names(coef(LVM1.red)))])
  
  method <- lava:::gaussian2_method.lvm
  suppressWarnings( 
    LVM2 <- estimate(m, data = d, control = list(iter.max = iter.max, start = start, method = method), estimator = "gaussian2")
  )
  # gaussian2_method.lvm
  suppressWarnings(
    LVM2.red <- estimate(m.red, data = d, control = list(iter.max = iter.max, start = start[coef(m.red)], method = method), estimator = "gaussian2")
  )
  
  expect_equal(coef(LVM2),coef(LVM2.red)[names(coef(LVM2))])
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

