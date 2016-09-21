path <- butils:::dir.gitHub()

library(lava)
library(testthat)
source(file.path(path,"lava.penalty/R/Lava_reducedMoment.R"))
source(file.path(path,"lava.penalty/R/Lava_reducedObject.R"))

m <- lvm()
m <- regression(m,y=c('y1','y2'),x='x'%++%1:20)
regression(m) <- y2~y1

# simul
set.seed(10)
d <- sim(m,30)

# normal model
e <- estimate(m,d)
logLik(m,data=d,p=coef(e))
score(m,data=d,p=coef(e),indiv=TRUE)
colSums(score(m,data=d,p=coef(e),indiv=TRUE))

# reduced model
mR <- lvm()
regression(mR) <- y2~y1  
mR <- regression(mR,y=c('y1','y2'),x='x'%++%1:20, reduce = TRUE)

# comparison of the moments
beta1 <- coef(e)
beta2 <- coef(e)+1

dcomp <- d
dcomp$LPy1 <- NA
dcomp$LPy2 <- NA

expect_equal(as.double(gaussianReduced_objective.lvm(mR, p = beta1[coef(mR)], data = dcomp)),
             as.double(logLik(m,data=d,p=beta1))
)
expect_equal(as.double(gaussianReduced_objective.lvm(mR, p = beta2[coef(mR)], data = dcomp)),
             as.double(logLik(m,data=d,p=beta2))
)
expect_equal(gaussianReduced_gradient.lvm(mR, p = beta1[coef(mR)], data = dcomp, indiv = TRUE),
             score(m,data=d,p=beta1,indiv=TRUE)[,coef(mR)]
)
expect_equal(gaussianReduced_gradient.lvm(mR, p = beta2[coef(mR)], data = dcomp, indiv = TRUE),
             score(m,data=d,p=beta2,indiv=TRUE)[,coef(mR)]
)
expect_equal(as.double(gaussianReduced_gradient.lvm(mR, p = beta1[coef(mR)], data = dcomp, indiv = FALSE)[coef(m)]),
             as.double(score(m,data=d,p=beta1,indiv=FALSE))
)
expect_equal(as.double(gaussianReduced_gradient.lvm(mR, p = beta2[coef(mR)], data = dcomp, indiv = FALSE)[coef(m)]),
             as.double(score(m,data=d,p=beta2,indiv=FALSE))
)

# estimation of the reduced model
res <- estimate(mR, data = d, control = list(start = beta1[coef(mR)], trace = 3),
                estimator = "gaussianReduced")
range(coef(res)-beta1[coef(mR)])

res0 <- estimate(m,d, control = list(iter.max = 0))
logLik(res0)

logLik(m,data=d,p=p[coef(m)])

start <- coef(res0)

####DO NOT WORK
resRed <- estimate(mR, data = d, control = list(trace = 3, start = start[coef(mR)], iter.max = 1),
                estimator = "gaussianReduced")
coef(resRed)

resNorm <- estimate(m, data = d, control = list(trace = 3, start = start, iter.max = 1))
(start-coef(resNorm))/derivSave

score(m,data=d,p=start)
score(m,data=d,p=p[coef(m)])

(start[coef(mR)]-p)/derivSave



m0 <- lvm(y1~1*LPy1,y2~y1+1*LPy2)

logL <- function(p) {
  ii <- match(coef(m0),names(p))
  p0 <- p[ii]
  iix <- setdiff(seq_along(p),ii)
  px <- p[iix]
  X1 <- as.matrix(d[,-c(1:2)])
  X2 <- as.matrix(d[,-c(1:2)])
  d0 <- transform(d,
                  lp1=as.vector(X2%*%px[1:20]),
                  lp2=as.vector(X1%*%px[21:40]))
  logLik(m0,data=d0,p=p0)
}


