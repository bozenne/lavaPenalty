### test-lassoRegression_proxGrad.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar  7 2017 (15:44) 
## Version: 
## last-updated: mar  9 2017 (14:21) 
##           By: Brice Ozenne
##     Update #: 28
#----------------------------------------------------------------------
## 
### Commentary:
##
## BEFORE RUNNING THE FILE:
## library(butils.base) ; package.source("lava.penalty") ;
## path <- file.path(butils.base::path_gitHub(),"lava.penalty","tests")
## source(file.path(path,"FCT.R"))
##
## CONTENT
## Test regression with only continuous variables: no penalty, at breakpoints, between breakpoints
##                 with partial penalization
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(lava.penalty)
library(penalized)
library(testthat)

#### > settings ####
test.tolerance <- 1e-4
test.scale <- NULL
lava.penalty.options(trace = FALSE, type = "lava")
lava.penalty.options(trace = FALSE, type = "proxGrad")

context("#### Estimate lasso regression with proximal gradient #### \n")

# {{{ simulation and LVM
# parametrisation
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( c(rep(0,2),1:3) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)

# simulation
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))

# estimation
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
# }}}
# {{{ no penalty
test_that("NR vs proxGrad with lasso - lambda=0", {
    plvm.model <- penalize(lvm.model)
    
    # likelihood
    eplvm.0 <- estimate(plvm.model, df.data, lambda1 = 0)
    expect_equal(object=coef(elvm.model),
                 expected=coef(eplvm.0),
                 tolerance=test.tolerance, scale=test.scale)    

    # Least square
    eplvm.0 <- estimate(plvm.model,
                        df.data,
                        lambda1 = 0,
                        constrain.lambda = TRUE)
               
    expect_equal(object=coef(elvm.model),
                 expected=coef(eplvm.0),
                 tolerance=test.tolerance,
                 scale=test.scale)    
})
# }}}

# {{{ proximal gradient algorithm at breakpoints

penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = FALSE)
seq_lambda1 <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seq_lambda1Sigma <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))
n.lambda <- length(seq_lambda1)

plvm.model <- penalize(lvm.model)

for(iter_l in 1:n.lambda){ # iter_l <- 5

    coef.expected <- coef2.penalized(penalized.PathL1[[iter_l]], name.response = "Y")[names(coef(elvm.model))]
    logLik.expected <- penalized.PathL1[[iter_l]]@loglik
    penalty.expected <- sum(penalized.PathL1[[iter_l]]@penalty)/penalized.PathL1[[iter_l]]@nuisance$sigma2
    iLambda1 <- seq_lambda1[iter_l]
    iLambda1Sigma <- seq_lambda1Sigma[iter_l]
    
    # normal model
    eplvm.fit_tempo1 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma,
                                 constrain.lambda = FALSE,
                                 control = list(constrain=FALSE))
    
    test_that("proxGrad with lasso", {
        expect_equal(object = coef(eplvm.fit_tempo1),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)    
    })
  
    # with constrains
    eplvm.fit_tempo2 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma,
                                 constrain.lambda = FALSE,
                                 control = list(constrain=TRUE))

  
    test_that("proxGrad with lasso (constrain.variance)", {
        expect_equal(object = coef(eplvm.fit_tempo2),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)  
    })

    # fixed sigma
    eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data,
                                 lambda1 = iLambda1,
                                 constrain.lambda = TRUE,
                                 control = list(constrain = TRUE))

    test_that("proxGrad with lasso (fixed sigma)", {
        expect_equal(object = coef(eplvm.fit_tempo3),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)  

        expect_equal(as.double(logLik(eplvm.fit_tempo3, addPenalty = FALSE)),
                     logLik.expected,
                     tolerance = test.tolerance, scale = 1)

        expect_equal(as.double(logLik(eplvm.fit_tempo3, addPenalty = TRUE)),
                     logLik.expected + penalty.expected,
                     tolerance = test.tolerance, scale = 1)
    })

    test_that("logLik for different parametrisations", {
        expect_equal(logLik(eplvm.fit_tempo1, addPenalty = FALSE),
                     logLik(eplvm.fit_tempo2, addPenalty = FALSE),
                     tolerance = test.tolerance, scale = 1)
        expect_equal(logLik(eplvm.fit_tempo1, addPenalty = TRUE),
                     logLik(eplvm.fit_tempo2, addPenalty = TRUE),
                     tolerance = test.tolerance, scale = 1)
    })
}

# }}}

# {{{ proximal gradient algorithm between breakpoints
seq2_lambda <- igraph::running.mean(seq_lambda, binwidth = 2)

for(iter_l in 1:length(seq2_lambda)){ # iter_l <- 1
    penalized.L1 <- penalized(Y ~  ., data = df.data, lambda1 = seq2_lambda[iter_l], trace = FALSE)

    lambda1_l <- penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2
    lambda1_lFIXED <- penalized.PathL1[[iter_l]]@lambda1

    # normal model
    eplvm.fit_tempo1 <- estimate(plvm.model, data = df.data, lambda1 = lambda1_l,
                                 constrain.lambda = FALSE,
                                 control = list(constrain=FALSE))
    
    test_that("NR vs proxGrad with lasso (between knots)", {
        expect_equal(object = coef(eplvm.fit_tempo1),
                     expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                     tolerance = test.tolerance, scale=1)    
    })
  
    # with constrains
    eplvm.fit_tempo2 <- estimate(plvm.model, data = df.data, lambda1 = lambda1_l,
                                 constrain.lambda = FALSE,
                                 control = list(constrain=TRUE))
  
    test_that("NR vs proxGrad with lasso (constrain.variance - between knots)", {
        expect_equal(object = coef(eplvm.fit_tempo2),
                     expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                     tolerance = test.tolerance, scale=1)  
    })

    # fixed sigma
    eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data,
                                 lambda1 = lambda1_lFIXED,
                                 constrain.lambda = TRUE, control = list(constrain = TRUE)
                                 )

    test_that("NR vs proxGrad with lasso (fix sigma - between knots)", {
        expect_equal(object = coef(eplvm.fit_tempo3),
                     expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                     tolerance = test.tolerance, scale=1)  
    })
}
# }}}

# {{{ partial penalization 
X.unpenalized <- cbind(1,as.matrix(df.data[,c("X2","X4")]))
X.penalized <- as.matrix(df.data[,c("X1","X3","X5")])

penalized.PathL1 <- penalized(response = df.data$Y,
                              unpenalized = X.unpenalized,
                              penalized = X.penalized,
                              steps = "Park", trace = FALSE
                              )

seq_lambda1 <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seq_lambda1Sigma <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))
n.lambda <- length(seq_lambda1)

plvm.model <- penalize(lvm.model, value = Y ~ X1 + X3 + X5)

for(iter_l in 1:n.lambda){ # iter_l <- 5

    coef.expected <- coef2.penalized(penalized.PathL1[[iter_l]], name.response = "Y")[names(coef(elvm.model))]
    logLik.expected <- penalized.PathL1[[iter_l]]@loglik
    penalty.expected <- sum(penalized.PathL1[[iter_l]]@penalty)/penalized.PathL1[[iter_l]]@nuisance$sigma2
    iLambda1 <- seq_lambda1[iter_l]
    iLambda1Sigma <- seq_lambda1Sigma[iter_l]
    
    # normal model
    eplvm.fit_tempo1 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma,
                                 constrain.lambda = FALSE,
                                 control = list(constrain=FALSE))
    
    test_that("proxGrad with lasso", {
        expect_equal(object = coef(eplvm.fit_tempo1),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)    
    })
  
    # with constrains
    eplvm.fit_tempo2 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma,
                                 constrain.lambda = FALSE,
                                 control = list(constrain=TRUE))

  
    test_that("proxGrad with lasso (constrain.variance)", {
        expect_equal(object = coef(eplvm.fit_tempo2),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)  
    })

    # fixed sigma
    eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data,
                                 lambda1 = iLambda1,
                                 constrain.lambda = TRUE,
                                 control = list(constrain = TRUE))

    test_that("proxGrad with lasso (fixed sigma)", {
        expect_equal(object = coef(eplvm.fit_tempo3),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)  

        expect_equal(as.double(logLik(eplvm.fit_tempo3, addPenalty = FALSE)),
                     logLik.expected,
                     tolerance = test.tolerance, scale = 1)

        expect_equal(as.double(logLik(eplvm.fit_tempo3, addPenalty = TRUE)),
                     logLik.expected + penalty.expected,
                     tolerance = test.tolerance, scale = 1)
    })

    test_that("logLik for different parametrisations", {
        expect_equal(logLik(eplvm.fit_tempo1, addPenalty = FALSE),
                     logLik(eplvm.fit_tempo2, addPenalty = FALSE),
                     tolerance = test.tolerance, scale = 1)
        expect_equal(logLik(eplvm.fit_tempo1, addPenalty = TRUE),
                     logLik(eplvm.fit_tempo2, addPenalty = TRUE),
                     tolerance = test.tolerance, scale = 1)
    })
}
# }}}

#----------------------------------------------------------------------
### test-lassoRegression_proxGrad.R ends here
