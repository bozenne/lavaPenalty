# library(ggplot2)
# library(lava)
library(lava.penalty)

EPSODE_options <- lava.options()$EPSODE
EPSODE_options$resolution_lambda1 <- c(1e-10,1)
lava.options(EPSODE = EPSODE_options)

####> linear regression ####
set.seed(10)
n <- 300
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:12), collapse = "+")))
mSim <- lvm(formula.lvm)
df.data <- sim(mSim,n)

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)

#### lasso penalty

## regularization path
path1F <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fixSigma = TRUE)

# backward
path1B <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE,
                   increasing = FALSE, fixSigma = TRUE)

getPath(path1B)

plot(path1B)
plot(path1B, type = "path")

## proxGrad
pfit <- estimate(plvm.model,  data = df.data,
                 lambda1 = getPath(path1B, names = "lambda1.abs")[6,1],
                 control = list(constrain = TRUE), fixSigma = TRUE)
pfit

#### ridge penalty
pfit <- estimate(plvm.model,  data = df.data, lambda2 = 10, control = list(constrain = TRUE))
pfit

#### elastic net penalty
pfit <- estimate(plvm.model,  data = df.data, lambda1 = 100, lambda2 = 10, control = list(constrain = TRUE))
pfit

#### group lasso penalty
plvm.model <- penalize(lvm.model, value = paste0("Y~X",1:4), group = 1)
plvm.model <- penalize(plvm.model, value = paste0("Y~X",5:9), group = 2)

pfit <- estimate(plvm.model,  data = df.data, lambda1 = 80, control = list(constrain = TRUE), fixSigma = TRUE)
pfit

#### nuclear norm penalty



####> Latent variable model ####

set.seed(10)
n <- 300
mSim <- lvm(list(Y1 ~ eta + X1, Y2 ~ eta, Y3 ~ eta + X2))
latent(mSim) <- ~eta
regression(mSim) <- eta~Z1
covariance(mSim) <- Y1~Y2

df.data <- sim(mSim,n)


#### partial lasso penalty
m <- mSim
regression(m) <- Y2~X1+X2
pm <- penalize(m, c("Y2~X1","Y2~X2","Y3~X2"))

## regularization path
path.lvm <- estimate(pm, df.data, regularizationPath = TRUE, stopParam = 1, increasing = FALSE)

pathLVM <- calcLambda(model = pm, data.fit = df.data, seq_lambda1 = seq(0,300, length.out = 5))
plot(pathLVM)
plot(pathLVM, type = "path")
# lava:::estimate.lvm(pm, scale(df.data))

## local estimation
proxGrad_options <- lava.options()$proxGrad
proxGrad_options$method <- "ISTA"
lava.options(proxGrad = proxGrad_options)

e <- estimate(pm, df.data, lambda1 = 0, control = list(constrain = TRUE))

proxGrad_options <- lava.options()$proxGrad
proxGrad_options$method <- "mFISTA"
lava.options(proxGrad = proxGrad_options)

e <- estimate(pm, df.data, lambda1 = 0, control = list(constrain = TRUE))


e <- estimate(pm, df.data, lambda1 = 2e1, control = list(constrain = TRUE))




