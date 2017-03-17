path <- file.path(butils::dir.gitHub(),"lava.penalty","tests")
path.res <- file.path(path, "Results/ElasticNet")

# library(butils)
# package.source("lava.penalty")
library(penalized)
library(data.table)
library(lava.penalty)
library(fields)
library(testthat)

context("#### Reg-NuclearNorm #### \n")

#### settings ####
set.seed(10)
n.obs <- 500
size <- "big"

#### 1- Simulation ###
xmax <- switch(size,
               small = 10,
               big = 64)
ymax <- switch(size,
               small = 10,
               big = 64)
center <- switch(size,
                 small = c(5,5),
                 big = c(32,32))
radius <- switch(size,
                 small = 3,
                 big = 10)
res <- simForm(n.obs, xmax, ymax, radius = radius)
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

## display
display <- FALSE
if(display){
fields:::image.plot(1:xmax, 1:ymax, res$X, 
                    breaks = seq(-1,1,0.05), col = tim.colors(40), 
                    xlab = "", ylab = "")}

#### 2- Model ####
names.param <- c("intercept",Znames,Xnames,"Y,Y")
beta <- setNames(c(0,gamma,betaI,1), names.param)
beta[Xnames] <-  0

#### lasso 
test.penaltyLasso <- setNames(names(beta)[-length(beta)] %in% c(Znames, Xnames), c(Znames, Xnames))
glmnet.fit <- glmnet:::glmnet(x = cbind(Z,X), y = Y, family = "gaussian", alpha = 0,
                              penalty.factor = test.penaltyLasso)

B.LS <- matrix(coef(glmnet.fit, s = 0)[-(1:(n.confounder+1))],
               nrow = xmax, ncol = ymax, byrow = TRUE, 
               breaks = seq(-1,1,0.05), col = tim.colors(40))
fields:::image.plot(B.LS,
                    breaks = seq(-1,1,0.05), col = tim.colors(40), 
                    xlab = "", ylab = "")

B.LS <- matrix(coef(glmnet.fit, s = 100)[-(1:(n.confounder+1))],
               nrow = xmax, ncol = ymax, byrow = TRUE)
fields:::image.plot(B.LS,
                    breaks = seq(-1,1,0.05), col = tim.colors(40), 
                    xlab = "", ylab = "")


loss.glmnet <- sapply(1:length(glmnet.fit$lambda), function(x){
  sum((glmnet.fit$beta[-(1:n.confounder),x]-betaI)^2)
})
plot(glmnet.fit$lambda,loss.glmnet)

#### nuclear norm
formula.image <- as.formula(paste0("Y~",paste(c(Znames,Xnames), collapse = "+")))

lvm.image <- lvm(as.formula(paste0("Y~",paste(Znames, collapse = "+"))))
plvm.image <- lvm.image
penalizeNuclear(plvm.image, coords = coords) <- as.formula(paste0("Y~",paste0(Xnames,collapse = "+")))

control.proxGrad <- lava.options()$proxGrad
control.proxGrad$method <- "ISTA"
lava.options(proxGrad = control.proxGrad)
elvm.Path <- estimate(plvm.image,  data = dt.data, lambdaN = 1e1, 
                      control = list(iter.max = 30, trace = 2))
B.LS <- matrix(attr(elvm.Path$opt$message,"par")[paste0("Y~",Xnames)],
               nrow = xmax, ncol = ymax, byrow = TRUE)
fields:::image.plot(B.LS)

elvm.Path <- estimate(plvm.image,  data = dt.data, lambdaN = 1e2,
                      control = list(iter.max = 30, trace = 2))
B.LS <- matrix(attr(elvm.Path$opt$message,"par")[paste0("Y~",Xnames)],
               nrow = xmax, ncol = ymax, byrow = TRUE)
fields:::image.plot(B.LS)

elvm.Path <- estimate(plvm.image,  data = dt.data, lambdaN = 5e2,
                      control = list(iter.max = 50, trace = 2))
B.LS <- matrix(attr(elvm.Path$opt$message,"par")[paste0("Y~",Xnames)],
               nrow = xmax, ncol = ymax, byrow = TRUE)
fields:::image.plot(B.LS)
