path <- butils:::dir.gitHub()

library(lava)
library(testthat)
# library(butils)
# package.source("lava", RorderDescription = FALSE)
source(file.path(path,"lava.penalty/R/Lava_reducedMoment.R"))
source(file.path(path,"lava.penalty/R/Lava_reducedObject.R"))

m <- lvm()
m <- regression(m,y=c('y1','y2'),x='x'%++%1:20)
regression(m) <- y2~y1

# simul
set.seed(10)
d <- sim(m,30)

# normal model
res0 <- lava:::estimate.lvm(m,d, control = list(iter.max = 0))
start <- coef(res0)


e <- estimate(m,d)
logLik(m,data=d,p=coef(e))
score(m,data=d,p=coef(e),indiv=TRUE)
colSums(score(m,data=d,p=coef(e),indiv=TRUE))

# reduced model
mR <- lvm()
regression(mR) <- y2~y1  
mR <- regression(mR,y=c('y1','y2'),x='x'%++%1:20, reduce = TRUE) # need to be at the end ...

mRR <- lvm()
regression(mRR) <- y2~y1  
mRR <- regression(mRR,y=c('y1','y2'),x='x'%++%1:20)
mRR <- regression(mRR,eta~Z1+Z2)
system.time(
  reduce(mRR)
)

#### comparison of the moments ####
indexRED <- match(coef(mR),coef(m))
beta1 <- coef(e)
beta2 <- start # coef(e)+1

dcomp <- d
dcomp$LPy1 <- NA
dcomp$LPy2 <- NA

for(iterB in 1:2){
  beta <- list(beta1,beta2)[[iterB]]
  ## test logLik
  expect_equal(as.double(gaussianReduced_logLik.lvm(mR, p = beta[indexRED], data = dcomp)),
               lava:::gaussian_logLik.lvm(object = m, data=d, p=beta)
  )
  
  ## test Objective
  # expect_equal(as.double(gaussianReduced_objective.lvm(mR, p = beta[coef(mR)], data = dcomp)),
  #              lava:::gaussian_objective.lvm(x = m, data=d, p=beta, n = e$data$n, S = e$S, mu = e$mu)
  # )
  
  ## test gradient
  expect_equal(unname(gaussianReduced_gradient.lvm(mR, p = beta[indexRED], data = dcomp, indiv = FALSE)),
               lava:::gaussian_gradient.lvm(x = m, data=d, p=beta, S = e$S, n = e$data$n, mu = e$mu)[indexRED]
  )
  
  ## test score
  expect_equal(unname(gaussianReduced_score.lvm(mR, p = beta[indexRED], data = dcomp, indiv = FALSE)),
               unname(lava:::gaussian_score.lvm(x = m, data=d, p=beta, S = e$S, n = e$data$n, mu = e$mu)[,indexRED,drop=FALSE])
  )
  expect_equal(gaussianReduced_score.lvm(mR, p = beta[indexRED], data = dcomp, indiv = TRUE),
               score(m,data=d,p=beta,indiv=TRUE)[,indexRED,drop=FALSE]
  )
  expect_equal(unname(gaussianReduced_score.lvm(mR, p = beta[indexRED], data = dcomp, indiv = FALSE)),
               unname(score(m,data=d,p=beta,indiv=FALSE)[,indexRED,drop=FALSE])
  )
  
  ## test hessian
  Hred <- gaussianReduced_hessian.lvm(x = mR, data=d, p=beta[indexRED], n = e$data$n, mu = e$mu, S = e$S, implementation = "test", type = "E")
  Hlava <- lava:::gaussian2_hessian.lvm(x = m, data=d, p=beta, n = e$data$n, mu = e$mu, S = e$S)
  expect_equal(unname(attr(Hred,"grad")),
               attr(Hlava,"grad")[indexRED]
  ) 
  attr(Hred,"grad") <- NULL # Hred[1:5,1:5]
  attr(Hlava,"grad") <- NULL # Hlava[1:10,1:10] # Hlava[indexRED,indexRED][1:10,1:10]
  expect_equal(unname(Hred), Hlava[indexRED,indexRED])

  Hred <- gaussianReduced_hessian.lvm(x = mR, data=d, p=beta[indexRED], n = e$data$n, mu = e$mu, S = e$S, implementation = "test", type = "num")
  Hlava <- lava:::gaussian1_hessian.lvm(x = m, data=d, p=beta, n = e$data$n, mu = e$mu, S = e$S)
  expect_equal(unname(Hred),Hlava[indexRED,indexRED], tolerance = 1e-3)
  cat("range of the difference in Hessian (numDeriv): ",paste(range(Hred-Hlava[indexRED,indexRED]),collapse = " "),"\n")
  
  ### do not work
  Hred <- gaussianReduced_hessian.lvm(x = mR, data=d, p=beta[indexRED], n = e$data$n, mu = e$mu, S = e$S, type = "information")
  Hlava <- lava:::gaussian_hessian.lvm(x = m, data=d, p=beta, n = e$data$n, mu = e$mu, S = e$S)
  # Hlava[indexRED[1:5],indexRED[1:5]]
  expect_equal(unname(Hred[1:5,1:5]),
               Hlava[indexRED,indexRED][1:5,1:5]
  )
  cat("range of the difference in Hessian (information): ",paste(range(Hred-Hlava[indexRED,indexRED]),collapse = " "),"\n")
  # fields:::image.plot(Hlava)
  # fields:::image.plot(Hred)
  # debug(lava:::information.lvm)
  
  ## test information
  # Ilava <- information(m,data=d,p=beta)
  
}  

#### estimation of the model ####
resGS <- estimate(m, data = d, control = list(trace = 3, start = coef(e)))
resRed <- estimate.lvm.reduced(mR, data = d, control = list(trace = 3, start = coef(e)[indexRED]),
                               estimator = "gaussianReduced")
expect_equal(coef(resRed),coef(resGS)[indexRED])

resLAVA <-  estimate(m, data = d, control = list(trace = 3, start = start))
coef(resLAVA)

resRed <- estimate.lvm.reduced(mR, data = d, control = list(trace = 3, start = start[indexRED], iter.max = 1000),
                estimator = "gaussianReduced")
coef(resRed)
#'log Lik.' -231.726 (df=7)

