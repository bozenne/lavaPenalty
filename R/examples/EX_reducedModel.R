if(debug){
  path.FCT <- butils:::dir.gitHub()
  
  # install_github("kkholst/lava", ref = "develop")
  library(testthat)
  library(butils)
  library(lava)
  # package.source("lava", RorderDescription = FALSE)
  package.source("lava.penalty", RorderDescription = FALSE)
  .onLoad() # gethook
}

#### regression ####

## simulation
m <- lvm()
m <- regression(m,y='y1',x='x'%++%1:2)
m <- regression(m,y='y1',x='z'%++%1:5)

set.seed(10)
d <- sim(m,100)

## reduced model 1
mR1 <- lvm()
mR1 <- regression(mR1,y='y1',x='x'%++%1:2)
mR1 <- regression(mR1,y='y1',x='z'%++%1:5, reduce = TRUE)

## reduced model 2
mR2 <- reduce(m)

## estimation
system.time(
  em1 <- estimate(m, d)
)
# system.time(
#   em2 <- estimate(m, d, estimator = "gaussian1")
# )
system.time(
  emR1 <- estimate(mR1, d)
)
coef(em1)-coef(emR1)[names(coef(em1))]


#### penalization ####
mp <- m
penalize(mp) <- grep("x", coef(m),value = TRUE)
rmp <- reduce(mp)

epm <- estimate(mp,d,lambda1 = 1e5, control = list (trace = 2, iter.max = 2000))
epmR <- estimate(rmp,d,lambda1 = 1e5, control = list (trace = 2, iter.max = 2000))


#### latent variable model ####
m <- lvm()
m <- regression(m,y=c('y1','y2','y3','y4'),x='eta')
m <- regression(m,y=c('y2','y3'),x='x'%++%1:5)
latent(m) <- ~eta
m <- regression(m,y=c('y1','y2'),x='z'%++%1:2)
covariance(m) <- y2~y1

# simul
set.seed(10)
d <- sim(m,100)

# reduced model
mR1 <- reduce(m, endo = c("y2"))
mR2 <- reduce(m)

## estimation
em <- estimate(m,d)
emR <- estimate(mR1, data = d, control = list(trace=2))

