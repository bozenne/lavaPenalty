path.FCT <- butils:::dir.gitHub()

# install_github("kkholst/lava", ref = "develop")
library(testthat)
library(butils)
library(lava)
# package.source("lava", RorderDescription = FALSE)
package.source("lava.penalty", RorderDescription = FALSE)
.onLoad() # gethook


#### regression ####
n.covar <- 3
n.z <- 2

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

#### H: numerical derivative
eNum <- estimate(m, d, 
                 estimator = "gaussian1",
                 control = list(method = "NR", trace = trace, iter.max = iter.max, constrain = TRUE))
coef(ebis)
eRNum <- estimate(mR,d, 
               estimator = "gaussian1",
               control = list(method = "NR", trace = trace, iter.max = iter.max, start = startE0[coef(mR)], constrain = TRUE))
eRNum <- estimate(mR,d, 
                  estimator = "gaussian1",
                  control = list(method = "NR", trace = trace, iter.max = iter.max, constrain = TRUE))

indexMatch <- match(names(coef(eRNum)),names(coef(eNum)))
expect_equal(coef(eRNum), coef(eNum)[indexMatch])
expect_equal(eRNum$vcov[1:length(indexMatch),1:length(indexMatch)], eNum$vcov[indexMatch,indexMatch])


#### H: S*t(S)
eS <- estimate(m, d, 
                 estimator = "gaussian2",
                 control = list(method = "NR", trace = trace, iter.max = iter.max))
eRS <- estimate(mR,d, 
               estimator = "gaussian2",
               control = list(method = "NR", trace = trace, iter.max = iter.max, start = startE0[coef(mR)]))

indexMatch <- match(names(coef(eRS)),names(coef(eS)))
expect_equal(coef(eRS), coef(eS)[indexMatch])
expect_equal(unname(eRS$vcov[1:length(indexMatch),1:length(indexMatch)]), eS$vcov[indexMatch,indexMatch])


#### H: information
eGS <- estimate(m, d, 
                estimator = "gaussian",
                control = list(method = "NR", trace = trace, iter.max = iter.max, constrain = TRUE))

expect_equal(coef(eGS), coef(eS))
expect_equal(coef(eGS), coef(eNum))


indexMatch <- match(names(coef(eR)),names(coef(e)))
expect_equal(coef(eR), coef(ebis)[indexMatch])
expect_equal(unname(eR$vcov[1:length(indexMatch),1:length(indexMatch)]), ebis$vcov[indexMatch,indexMatch])



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

## without second moment
resGS <- estimate(m, data = d, control = list(trace = 3, start = startE, method = "NR"))
resRed <- estimate.lvm.reduced(mR, data = d, control = list(trace = 3, start = startRE, method = "NR"))

resGS <- estimate(m, data = d, control = list(trace = 3, start = start0))
resRed <- estimate.lvm.reduced(mR, data = d, control = list(trace = 3, start = startR0))
expect_equal(coef(resRed),coef(resGS)[indexRED])

resGS <- estimate(m, data = d, control = list(trace = 3, start = start0), method = "gaussian")


#### link with penalization ####
n.covar <- 10
n.z <- 2 

m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:n.covar)
m <- regression(m,y='y1',x='z'%++%1:n.z)

# simul
set.seed(10)
d <- sim(m,100)

# normal model
mp <- m
penalize(mp) <- grep("x", coef(m),value = TRUE)
mp <- reduce(mp)

e <- estimate.plvm(mp,d,lambda1 = 2, lambda2 = 10)

lv1 <- lvm(Y ~ X1 + X2)
lv2 <- lv1
parameter(lv2, remove = FALSE) <- "Y~X3"
coef(lv2)
coef(lv1)

res <- rmLink(lv2,"Y~X1")
coef(res)

df <- sim(lv1, n = 1e3)
res1 <- estimate(lv1, data = df)
res2 <- estimate.lvm(lv2, data = df)

coef(res1)
coef(res2)


# debug(`parameter<-`)
