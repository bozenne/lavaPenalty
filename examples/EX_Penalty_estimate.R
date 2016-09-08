# library(ggplot2)
# library(lava)

####> linear regression ####
set.seed(10)
n <- 300
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:12), collapse = "+")))
mSim <- lvm(formula.lvm)
categorical(mSim,labels=c("A","B","C")) <- "X1"
categorical(mSim,labels=c("A","B","C")) <- "X11"
df.data <- sim(mSim,n)

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)

## no penalty
fit <- estimate(lvm.model,  data = df.data)
fit

pfit <- estimate(plvm.model,  data = df.data, lambda1 = 0, control = list(constrain = TRUE))
pfit

#### lasso penalty
pfit <- estimate(plvm.model,  data = df.data, lambda1 = 116, control = list(constrain = TRUE))
pfit

## regularization path
path1F <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fixSigma = TRUE)
getPath(path1F, order = "lambda1.abs")

path1B <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fixSigma = TRUE, 
                   increasing = FALSE, resolution_lambda1 = c(1,1e-3))
path1B

getPath(path1B, order = "lambda1.abs")

#### ridge penalty


#### elastic net penalty

#### group lasso penalty

#### nuclear norm penalty


plot(path1B, type = "path", lambda = "lambda1.abs")
plot(path1B,  lambda = "lambda1.abs")
