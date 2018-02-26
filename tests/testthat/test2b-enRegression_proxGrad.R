### test-ridgeRegression_proxGrad.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar 15 2017 (17:39) 
## Version: 
## last-updated: mar 16 2017 (09:38) 
##           By: Brice Ozenne
##     Update #: 12
#----------------------------------------------------------------------
## 
### Commentary:
##
## BEFORE RUNNING THE FILE:
## library(butils.base) ; package.source("lavaPenalty") ;
## path <- file.path(butils.base::path_gitHub(),"lavaPenalty","tests")
## source(file.path(path,"FCT.R"))
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(lavaPenalty)
library(penalized)
library(testthat)

#### > settings ####
test.tolerance <- 1e-4
test.scale <- NULL
lavaPenalty.options(trace = FALSE, type = "lava")
lavaPenalty.options(trace = FALSE, type = "proxGrad")
lambda2 <- 2

context("#### Estimate elasticNet regression with proximal gradient #### \n")

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

# {{{ proximal gradient algorithm at breakpoints

penalized.PathL1 <- penalized(Y ~  ., data = df.data, lambda2 = lambda2, steps = "Park", trace = FALSE)
seq_lambda1 <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seq_lambda1Sigma <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))
seq_lambda2Sigma <- unlist(lapply(penalized.PathL1, function(x){lambda2/x@nuisance$sigma2}))
n.lambda <- length(seq_lambda1)

plvm.model <- penalize(lvm.model)

for(iter_l in 1:n.lambda){ # iter_l <- 4

    coef.expected <- coef2.penalized(penalized.PathL1[[iter_l]], name.response = "Y")[names(coef(elvm.model))]
    iLambda1 <- seq_lambda1[iter_l]
    iLambda1Sigma <- seq_lambda1Sigma[iter_l]
    iLambda2 <- lambda2
    iLambda2Sigma <- seq_lambda2Sigma[iter_l]

    logLik.expected <- penalized.PathL1[[iter_l]]@loglik
    penalty.expected1 <- sum(abs(penalized.PathL1[[iter_l]]@penalized))
    penalty.expected2 <- sum(penalized.PathL1[[iter_l]]@penalized^2)
    penalty.expected <- iLambda1Sigma * penalty.expected1 + iLambda2Sigma/2*penalty.expected2

    # normal model
    eplvm.fit_tempo1 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma, lambda2 = iLambda2Sigma,
                                 equivariance = FALSE,
                                 control = list(constrain=FALSE))

    test_that("proxGrad with lasso", {
        expect_equal(object = coef(eplvm.fit_tempo1),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)    
    })
  
    # with constrains
    eplvm.fit_tempo2 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma, lambda2 = iLambda2Sigma,
                                 equivariance = FALSE,
                                 control = list(constrain=TRUE))

  
    test_that("proxGrad with lasso (constrain.variance)", {
        expect_equal(object = coef(eplvm.fit_tempo2),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)  
    })

    # fixed sigma
    eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data,
                                 lambda1 = iLambda1, lambda2 = lambda2,
                                 equivariance = TRUE,#TRUE,
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
}

# }}}


# {{{ proximal gradient algorithm between breakpoints
seq2_lambda <- igraph::running.mean(seq_lambda1, binwidth = 2)

for(iter_l in 1:length(seq2_lambda)){ # iter_l <- 3
    penalized.L1 <- penalized(Y ~  ., data = df.data,
                              lambda1 = seq2_lambda[iter_l], lambda2 = lambda2,
                              trace = FALSE)

    coef.expected <- coef2.penalized(penalized.L1, name.response = "Y")[names(coef(elvm.model))]
    iLambda1 <- penalized.L1@lambda1
    iLambda1Sigma <- penalized.L1@lambda1/penalized.L1@nuisance$sigma2
    iLambda2 <- lambda2
    iLambda2Sigma <- lambda2/penalized.L1@nuisance$sigma2

    logLik.expected <- penalized.L1@loglik
    penalty.expected1 <- sum(abs(penalized.L1@penalized))
    penalty.expected2 <- sum(penalized.L1@penalized^2)
    penalty.expected <- iLambda1Sigma * penalty.expected1 + iLambda2Sigma/2*penalty.expected2
    
    # normal model
    eplvm.fit_tempo1 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma, lambda2 = iLambda2Sigma,
                                 equivariance = FALSE,
                                 control = list(constrain=FALSE))
    
    test_that("NR vs proxGrad with lasso (between knots)", {
        expect_equal(object = coef(eplvm.fit_tempo1),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)    
    })
  
    # with constrains
    eplvm.fit_tempo2 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma, lambda2 = iLambda2Sigma,
                                 equivariance = FALSE, control = list(constrain=TRUE))
    
    ## test_that("NR vs proxGrad with lasso (constrain.variance - between knots)", {
    ##     expect_equal(object = coef(eplvm.fit_tempo2),
    ##                  expected = coef.expected,
    ##                  tolerance = test.tolerance, scale=1)  
    ## })

    # fixed sigma
    eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data,
                                 lambda1 = iLambda1, lambda2 = iLambda2,
                                 equivariance = TRUE, control = list(constrain = TRUE)
                                 )
    
    test_that("NR vs proxGrad with lasso (fix sigma - between knots)", {
        expect_equal(object = coef(eplvm.fit_tempo3),
                     expected = coef.expected,
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
                              lambda2 = lambda2,
                              steps = "Park", trace = FALSE
                              )

seq_lambda1 <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seq_lambda1Sigma <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))
seq_lambda2Sigma <- unlist(lapply(penalized.PathL1, function(x){lambda2/x@nuisance$sigma2}))
n.lambda <- length(seq_lambda1)

plvm.model <- penalize(lvm.model, value = Y ~ X1 + X3 + X5)

for(iter_l in 1:n.lambda){ # iter_l <- 4

    coef.expected <- coef2.penalized(penalized.PathL1[[iter_l]], name.response = "Y")[names(coef(elvm.model))]
    iLambda1 <- seq_lambda1[iter_l]
    iLambda1Sigma <- seq_lambda1Sigma[iter_l]
    iLambda2 <- lambda2
    iLambda2Sigma <- seq_lambda2Sigma[iter_l]

    logLik.expected <- penalized.PathL1[[iter_l]]@loglik
    penalty.expected1 <- sum(abs(penalized.PathL1[[iter_l]]@penalized))
    penalty.expected2 <- sum(penalized.PathL1[[iter_l]]@penalized^2)
    penalty.expected <- iLambda1Sigma * penalty.expected1 + iLambda2Sigma/2*penalty.expected2
    
    # normal model
    eplvm.fit_tempo1 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma, lambda2 = iLambda2Sigma,
                                 equivariance = FALSE,
                                 control = list(constrain=FALSE))
    
    test_that("proxGrad with lasso", {
        expect_equal(object = coef(eplvm.fit_tempo1),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)    
    })
  
    # with constrains
    eplvm.fit_tempo2 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma, lambda2 = iLambda2Sigma,
                                 equivariance = FALSE,
                                 control = list(constrain=TRUE))
  
    test_that("proxGrad with lasso (constrain.variance)", {
        expect_equal(object = coef(eplvm.fit_tempo2),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)  
    })

    # fixed sigma
    eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data,
                                 lambda1 = iLambda1, lambda2 = iLambda2,
                                 equivariance = TRUE,
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
### test-ridgeRegression_proxGrad.R ends here
