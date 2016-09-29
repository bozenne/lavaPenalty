path.FCT <- butils:::dir.gitHub()

library(lava) # install_github("kkholst/lava", ref = "develop")
library(testthat)
library(butils)
package.source("lava", RorderDescription = FALSE)
source(file.path(path.FCT,"lava.penalty/R/Lava_reducedMoment.R"))
source(file.path(path.FCT,"lava.penalty/R/Lava_reducedObject.R"))
source(file.path(path.FCT,"lava.penalty/R/Lava_tempo.R"))



#### regression ####
n.covar <- 10
n.z <- 1

m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:n.covar)
m <- regression(m,y='y1',x='z'%++%1:n.z)

# simul
set.seed(10)
d <- sim(m,100)

# normal model
e0 <- lava:::estimate.lvm(m,d, control = list(iter.max = 0))
startE0 <- coef(e0)

e <- lava:::estimate.lvm(m,d)

# reduced model
mR <- lvm()
mR <- regression(mR,y='y1',x='x'%++%1:n.covar, reduce = TRUE)
mR <- regression(mR,y='y1',x='z'%++%1:n.z)

test <- checkMoment(mR, e, param = startE0)
test <- checkMoment(mR, e)

test$Hlava-test$Hred_num

system.time(
eR <- estimate(mR,d, 
               control = list(method = "NR", trace = 0, iter.max = 50, start = startE0[coef(mR)]))
)
system.time(
  e <- estimate.lvm(m,d, control = list(method = "NR", trace = 0, iter.max = 50), estimator = "gaussian1")
) 
system.time(
  eGS <- estimate.lvm(m,d)
)


eR <- estimate(mR,d, 
               control = list(method = "NR", trace = 0, iter.max = 50, start = startE0[coef(mR)]),
               estimator = "gaussianReduced1"
               )


range(coef(eR) - coef(e)[names(coef(eR))])
range(coef(eR) - coef(eGS)[names(coef(eR))])
range(vcov(eR) - vcov(e))
range(vcov(eR) - vcov(eGS))


 
#### latent variable model ####

m <- lvm()
m <- regression(m,y=c('y2','y3'),x='x'%++%1:5)
m <- regression(m,y=c('y1','y2','y3'),x='eta')
latent(m) <- ~eta
m <- regression(m,y=c('y1','y2'),x='z'%++%1:2)
covariance(m) <- y2~y1

# simul
set.seed(10)
d <- sim(m,100)

# normal model
res0 <- lava:::estimate.lvm(m,d, control = list(iter.max = 0))
start0 <- coef(res0)

e <- lava:::estimate.lvm(m,d)
startE <- coef(e)

logLik(e,data=d,p=coef(e))
colSums(score(e,data=d,p=coef(e),indiv=TRUE))

# reduced model
mR <- lvm()
mR <- regression(mR,y=c('y3','y2'),x='x'%++%1:5, reduce = TRUE)
mR <- regression(mR,y=c('y1','y2','y3'),x='eta')
latent(mR) <- ~eta
mR <- regression(mR,y=c('y1'),x='z'%++%1:2)
covariance(mR) <- y2~y1  
startR0 <- start0[match(coef(mR),coef(m))]
startRE <- startE[match(coef(mR),coef(m))]

system.time(
  mRR <- reduce(m)
)

#### estimation of the model ####
source(file.path(path.FCT,"lava.penalty/R/Lava_reducedMoment.R"))
source(file.path(path.FCT,"lava.penalty/R/Lava_reducedObject.R"))
source(file.path(path.FCT,"lava.penalty/R/Lava_tempo.R"))

## without second moment
resGS <- estimate(m, data = d, control = list(trace = 3, start = startE, method = "NR"))
resRed <- estimate.lvm.reduced(mR, data = d, control = list(trace = 3, start = startRE, method = "NR"))

resGS <- estimate(m, data = d, control = list(trace = 3, start = start0))
resRed <- estimate.lvm.reduced(mR, data = d, control = list(trace = 3, start = startR0))
expect_equal(coef(resRed),coef(resGS)[indexRED])

resGS <- estimate(m, data = d, control = list(trace = 3, start = start0), method = "gaussian")


#### comparison of the moments ####
# length(coef(mR)) - length(coef(m))
resR <- estimate(mR, data = d, control = list(trace = 3, iter.max = 0))

setdiff(coef(mR2), coef(mR))


beta1 <- startE
beta2 <- start0 # coef(e)+1

for(iterB in 1:2){
  beta <- list(beta1,beta2)[[iterB]]
  indexRED <- na.omit(match(coef(mR),names(beta)))
  ## test logLik
  expect_equal(gaussianReduced_logLik.lvm(mR, p = beta[indexRED], data = d),
               lava:::gaussian_logLik.lvm(object = e, data=d, p=beta) # normal_objective.lvm(e, data=d, p=beta)
  )
  
  ## test Objective
  # expect_equal(as.double(gaussianReduced_objective.lvm(mR, p = beta[coef(mR)], data = dcomp)),
  #              lava:::gaussian_objective.lvm(x = m, data=d, p=beta, n = e$data$n, S = e$S, mu = e$mu)
  # )
  
  ## test gradient
  expect_equal(unname(gaussianReduced_gradient.lvm(mR, p = beta[indexRED], data = d, indiv = FALSE)),
               lava:::gaussian_gradient.lvm(x = m, data=d, p=beta, S = e$S, n = e$data$n, mu = e$mu)[indexRED]
  )
  
  ## test score
  expect_equal(unname(gaussianReduced_score.lvm(mR, p = beta[indexRED], data = d, indiv = FALSE)),
               unname(lava:::gaussian_score.lvm(x = m, data=d, p=beta, S = e$S, n = e$data$n, mu = e$mu)[,indexRED,drop=FALSE])
  )
  expect_equal(gaussianReduced_score.lvm(mR, p = beta[indexRED], data = d, indiv = TRUE),
               score(m,data=d,p=beta,indiv=TRUE)[,indexRED,drop=FALSE]
  )
  expect_equal(unname(gaussianReduced_score.lvm(mR, p = beta[indexRED], data = d, indiv = FALSE)),
               unname(score(m,data=d,p=beta,indiv=FALSE)[,indexRED,drop=FALSE])
  )
  
  ## test hessian
  Hred <- gaussianReduced_hessian.lvm(x = mR, data=d, p=beta[indexRED], n = e$data$n,  type = "E")
  Hlava <- lava:::gaussian2_hessian.lvm(x = m, data=d, p=beta, n = e$data$n, mu = e$mu, S = e$S)
  expect_equal(attr(Hred,"grad"),
               attr(Hlava,"grad")[indexRED]
  ) 
  attr(Hred,"grad") <- NULL
  attr(Hlava,"grad") <- NULL
  expect_equal(Hred, Hlava[indexRED,indexRED])

  Hred <- gaussianReduced_hessian.lvm(x = mR, data=d, p=beta[indexRED], type = "num")
  Hlava <- lava:::gaussian1_hessian.lvm(x = m,  p=beta, n = e$data$n, mu = e$mu, S = e$S)
  expect_equal(Hred,Hlava[indexRED,indexRED])
  
  ### do not work
  Hred <- gaussianReduced_hessian.lvm(x = mR, data=d, p=beta[indexRED], type = "information")
  Hlava <- lava:::gaussian_hessian.lvm(x = m, data=d, p=beta, n = e$data$n, mu = e$mu, S = e$S)
  Hred-Hlava[indexRED,indexRED]
  # Hlava[indexRED[1:5],indexRED[1:5]]

  # expect_equal(unname(Hred[1:5,1:5]),
  #              Hlava[indexRED,indexRED][1:5,1:5]
  # )
  cat("range of the difference in Hessian (information): ",paste(range(Hred-Hlava[indexRED,indexRED]),collapse = " "),"\n")
  # fields:::image.plot(Hlava)
  # fields:::image.plot(Hred)
  
  ## test information
  # Ilava <- information(m,data=d,p=beta)
  
}  


