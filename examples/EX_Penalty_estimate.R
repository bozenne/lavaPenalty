# library(ggplot2)
# library(lava)
library(lava.penalty)
library(ggplot2)

EPSODE_options <- lava.options()$EPSODE
EPSODE_options$resolution_lambda1 <- c(1e-10,1)
lava.options(EPSODE = EPSODE_options)
options.proxGrad <- lava.options()$proxGrad
options.proxGrad$export.iter <- TRUE
lava.options(proxGrad = options.proxGrad)

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

#### algorithm ISTA
lambda1 <- 10

pfit_ISTA <- estimate(plvm.model,  data = df.data, method.proxGrad = "ISTA",
                      lambda1 = lambda1, iter.max = 100,
                      control = list(constrain = TRUE), fixSigma = TRUE)

pfit_FISTA_Vand <- estimate(plvm.model,  data = df.data, method.proxGrad = "FISTA_Beck",
                            lambda1 = lambda1, iter.max = 100,
                            control = list(constrain = TRUE), fixSigma = TRUE)

pfit_FISTA_Beck <- estimate(plvm.model,  data = df.data, method.proxGrad = "FISTA_Vand",
                            lambda1 = lambda1, iter.max = 100,
                            control = list(constrain = TRUE), fixSigma = TRUE)

pfit_mFISTA_Vand <- estimate(plvm.model,  data = df.data, method.proxGrad = "mFISTA_Vand",
                             lambda1 = lambda1, iter.max = 100,
                             control = list(constrain = TRUE), fixSigma = TRUE)


df.cv <- rbind(data.frame(pfit_ISTA$opt$details.cv, method = "ISTA"),
               data.frame(pfit_FISTA_Vand$opt$details.cv, method = "FISTA_Beck"),
               data.frame(pfit_FISTA_Beck$opt$details.cv, method = "FISTA_Vand"),
               data.frame(pfit_mFISTA_Vand$opt$details.cv, method = "mFISTA_Vand"))
gg <- ggplot(df.cv, aes(x = iteration, y = obj, group = method, color = method)) 
gg <- gg + geom_line() + geom_point()
gg

 gg + coord_cartesian(ylim = c(1490,1510), xlim = c(0,100))
gg + coord_cartesian(ylim = c(1500,1510), xlim = c(0,100))

#### nuclear norm penalty
n.obs <- 100
res <- simForm(n.obs, xmax = 25, ymax = 25, radius = 5)
coords <- res$coords
n.coord <- nrow(coords)
betaI <- as.vector(res$X)
X <- matrix(rnorm(n.obs*n.coord), nrow = n.obs, ncol = n.coord)
Xnames <- paste0("X",1:n.coord)

n.confounder <- 5
gamma <- rep(1, n.confounder)
Z <- matrix(rnorm(n.obs*n.confounder), nrow = n.obs, ncol = n.confounder)
Znames <- paste0("Z",1:n.confounder)

Y <- Z %*% gamma + X %*% betaI + rnorm(n.obs)
formula.lvm <- as.formula( paste0("Y~", paste0("X",1:n.coord, collapse = "+") ) )

dt.data <- data.table(Y=Y,data.frame(Z),data.frame(X))
names(dt.data) <- c("Y",Znames,Xnames)
dt.data[, (names(dt.data)) := lapply(.SD,scale), .SDcols = names(dt.data)]

##
lvm.image <- lvm(as.formula(paste0("Y~",paste(Znames, collapse = "+"))))
plvm.image <- lvm.image
penalizeNuclear(plvm.image, coords = coords) <- as.formula(paste0("Y~",paste0(Xnames,collapse = "+")))

for(lambda in c(1e0,1e1,1e2,2e2)){
  elvm.Path <- estimate(plvm.image,  data = dt.data, lambdaN = lambda, method.proxGrad = "FISTA_Vand",
                        control = list(iter.max = 100, constrain = TRUE))
  B.LS <- matrix(attr(elvm.Path$opt$message,"par")[paste0("Y~",Xnames)],
                 nrow = NROW(res$X), ncol = NCOL(res$X), byrow = TRUE)
  fields:::image.plot(B.LS)
  
  plot(elvm.Path$opt$details.cv[,"obj"])
  
}


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

#### algorithm ISTA
lambda1 <- 10

e_ISTA <- estimate(pm,  data = df.data, method.proxGrad = "ISTA",
                      lambda1 = lambda1,
                      control = list(constrain = TRUE))

e_FISTA_Vand <- estimate(pm,  data = df.data, method.proxGrad = "FISTA_Beck",
                            lambda1 = lambda1,
                            control = list(constrain = TRUE))

e_FISTA_Beck <- estimate(pm,  data = df.data, method.proxGrad = "FISTA_Vand",
                            lambda1 = lambda1,
                            control = list(constrain = TRUE))

e_mFISTA_Vand <- estimate(pm,  data = df.data, method.proxGrad = "mFISTA_Vand",
                             lambda1 = lambda1,
                             control = list(constrain = TRUE))


df.cv <- rbind(data.frame(e_ISTA$opt$details.cv, method = "ISTA"),
               data.frame(e_FISTA_Vand$opt$details.cv, method = "FISTA_Beck"),
               data.frame(e_FISTA_Beck$opt$details.cv, method = "FISTA_Vand"),
               data.frame(e_mFISTA_Vand$opt$details.cv, method = "mFISTA_Vand"))
gg <- ggplot(df.cv, aes(x = iteration, y = obj, group = method, color = method)) 
gg <- gg + geom_line() + geom_point() +
gg
