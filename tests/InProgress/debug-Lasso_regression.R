library(penalized)
library(testthat)
library(lava.penalty)
# butils::package.source("lava.penalty")

#### simulation ####
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( c(rep(0,2),1:3) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))

#### models  ####
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE, lambda2 = 0)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

# derivatives
lvGaussianO <- function(x){lvGaussian(x, df.data[[1]], as.matrix(df.data[,-1]))}
scoreGaussianO <- function(x){scoreGaussian(x, df.data[[1]], as.matrix(df.data[,-1]))}
hessianGaussianO <- function(x){hessianGaussian(x, df.data[[1]], as.matrix(df.data[,-1]))}

#### test ####

## at fix Sigma - forward
system.time(
  P1 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, fixSigma = TRUE,
                 regularizationPath = TRUE)
)

system.time(
  P2 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, fixSigma = TRUE,
                 objective = lvGaussianO, gradient = scoreGaussianO, hessian = hessianGaussianO, 
                 regularizationPath = TRUE)
)
# identical
expect_equal(P2$regPath, getPath(P1), tolerance = 1e-5)

## with free Sigma - forward
P3 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL,
               regularizationPath = TRUE)
getPath(P3)


P4 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, exportAllPath = FALSE,
               objective = lvGaussianO, gradient = scoreGaussianO, hessian = hessianGaussianO, 
               resolution_lambda1 = c(1e-10,1), # tries every lambda
               regularizationPath = TRUE)
getPath(P4$regPath)
# getPath(P1) - getPath(P4$regPath)

P5 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, 
               objective = lvGaussianO, gradient = scoreGaussianO, hessian = hessianGaussianO, 
               resolution_lambda1 = c(1e-10,0.1), # tries every lambda
               regularizationPath = TRUE)
getPath(P5$regPath)
# getPath(P1) - getPath(P5$regPath)

P6 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, exportAllPath = TRUE, 
               objective = lvGaussianO, gradient = scoreGaussianO, hessian = hessianGaussianO, 
               resolution_lambda1 = c(1e-10,0.01), # tries every lambda
               regularizationPath = TRUE)
getPath(P6$regPath)
# getPath(P1) - getPath(P6$regPath)

P7 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, 
               objective = lvGaussianO, gradient = scoreGaussianO, hessian = hessianGaussianO, 
               resolution_lambda1 = c(1e-2,1e-3), # adaptive mesh
               regularizationPath = TRUE)
getPath(P7$regPath)
# getPath(P1) - getPath(P7$regPath)

#### debug P6
n.points <- NROW(P6$regPath$path)
plot(P6$regPath, lambda = "lambda1", getCoef = "all",
     row = seq(1, n.points, length = 1000)) + coord_cartesian(ylim =c(-0.5,1.25), xlim = c(0,375))

plot(P6$regPath, lambda = "lambda1", getCoef = "all",
     row = seq(1, n.points, length = 1000)) + coord_cartesian(ylim = c(-0.5,1.25), xlim = c(340,350))
