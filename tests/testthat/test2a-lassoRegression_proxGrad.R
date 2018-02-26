### test-lassoRegression_proxGrad.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar  7 2017 (15:44) 
## Version: 
## last-updated: apr 20 2017 (15:55) 
##           By: Brice Ozenne
##     Update #: 44
#----------------------------------------------------------------------
## 
### Commentary:
##
## BEFORE RUNNING THE FILE:
## library(butils.base) ; package.source("lavaPenalty") ;
## path <- file.path(butils.base::path_gitHub(),"lavaPenalty","tests")
## source(file.path(path,"FCT.R"))
##
## CONTENT
## Test regression with only continuous variables: no penalty, at breakpoints, between breakpoints
##                 with partial penalization
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
    plvm.model <- penalize(lvm.model, Vridge = FALSE)
    
    # likelihood
    eplvm.0 <- estimate(plvm.model, df.data, lambda1 = 0)

    expect_equal(object=coef(elvm.model),
                 expected=coef(eplvm.0),
                 tolerance=test.tolerance, scale=test.scale)    

    # Least square
    eplvm.0 <- estimate(plvm.model,
                        df.data,
                        lambda1 = 0,
                        equivariance = TRUE)
               
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

plvm.model <- penalize(lvm.model, Vridge = FALSE)

for(iter_l in 1:n.lambda){ # iter_l <- 5

    coef.expected <- coef2.penalized(penalized.PathL1[[iter_l]], name.response = "Y")[names(coef(elvm.model))]
    logLik.expected <- penalized.PathL1[[iter_l]]@loglik
    penalty.expected <- sum(penalized.PathL1[[iter_l]]@penalty)/penalized.PathL1[[iter_l]]@nuisance$sigma2
    iLambda1 <- seq_lambda1[iter_l]
    iLambda1Sigma <- seq_lambda1Sigma[iter_l]
    
    # normal model
    eplvm.fit_tempo1 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma,
                                 equivariance = FALSE,
                                 control = list(constrain=FALSE))
    
    test_that("proxGrad with lasso", {
        expect_equal(object = coef(eplvm.fit_tempo1),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)    
    })
  
    # with constrains
    eplvm.fit_tempo2 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma,
                                 equivariance = FALSE,
                                 control = list(constrain=TRUE))

  
    test_that("proxGrad with lasso (constrain.variance)", {
        expect_equal(object = coef(eplvm.fit_tempo2),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)  
    })

    # fixed sigma
    eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data,
                                 lambda1 = iLambda1,
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
                              lambda1 = seq2_lambda[iter_l],
                              trace = FALSE)

    coef.expected <- coef2.penalized(penalized.L1, name.response = "Y")[names(coef(elvm.model))]
    iLambda1 <- penalized.L1@lambda1
    iLambda1Sigma <- penalized.L1@lambda1/penalized.L1@nuisance$sigma2

    logLik.expected <- penalized.L1@loglik
    penalty.expected <- iLambda1Sigma * sum(abs(penalized.L1@penalized))
    
    # normal model
    eplvm.fit_tempo1 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma, 
                                 equivariance = FALSE,
                                 control = list(constrain=FALSE))
    
    ## test_that("NR vs proxGrad with lasso (between knots)", {
    ##     expect_equal(object = coef(eplvm.fit_tempo1),
    ##                  expected = coef.expected,
    ##                  tolerance = test.tolerance, scale=1)    
    ## })
  
    # with constrains
    eplvm.fit_tempo2 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma, 
                                 equivariance = FALSE, control = list(constrain=TRUE))
    
    ## test_that("NR vs proxGrad with lasso (constrain.variance - between knots)", {
    ##     expect_equal(object = coef(eplvm.fit_tempo2),
    ##                  expected = coef.expected,
    ##                  tolerance = test.tolerance, scale=1)  
    ## })

    # fixed sigma
    eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data,
                                 lambda1 = iLambda1,
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
                              steps = "Park", trace = FALSE
                              )

seq_lambda1 <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seq_lambda1Sigma <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))
n.lambda <- length(seq_lambda1)

plvm.model <- penalize(lvm.model, value = Y ~ X1 + X3 + X5, Vridge = FALSE)

for(iter_l in 1:n.lambda){ # iter_l <- 5

    coef.expected <- coef2.penalized(penalized.PathL1[[iter_l]], name.response = "Y")[names(coef(elvm.model))]
    logLik.expected <- penalized.PathL1[[iter_l]]@loglik
    penalty.expected <- sum(penalized.PathL1[[iter_l]]@penalty)/penalized.PathL1[[iter_l]]@nuisance$sigma2
    iLambda1 <- seq_lambda1[iter_l]
    iLambda1Sigma <- seq_lambda1Sigma[iter_l]
    
    # normal model
    eplvm.fit_tempo1 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma,
                                 equivariance = FALSE,
                                 control = list(constrain=FALSE))
    
    test_that("proxGrad with lasso", {
        expect_equal(object = coef(eplvm.fit_tempo1),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)    
    })
  
    # with constrains
    eplvm.fit_tempo2 <- estimate(plvm.model, data = df.data,
                                 lambda1 = iLambda1Sigma,
                                 equivariance = FALSE,
                                 control = list(constrain=TRUE))

  
    test_that("proxGrad with lasso (constrain.variance)", {
        expect_equal(object = coef(eplvm.fit_tempo2),
                     expected = coef.expected,
                     tolerance = test.tolerance, scale=1)  
    })

    # fixed sigma
    eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data,
                                 lambda1 = iLambda1,
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


# {{{ multigroup

## simulate interaction
set.seed(10)
mSim <- lvm(Y ~ X1 + gender + group + Interaction)
distribution(mSim, ~gender) <- binomial.lvm()
distribution(mSim, ~group) <- binomial.lvm(size = 3)
constrain(mSim, Interaction ~ gender + group) <- function(x){x[,1]*x[,2]}
d <- sim(mSim, 1e2)
#d$gender <- factor(d$gender, labels = letters[1:2])
#d$group <- factor(d$group)

m <- lvm(Y ~ X1 + group)

test_that("same coef", {
    e <- estimate(list(m,m), split(d,d$gender))
    pm <- penalize(m, Vridge = FALSE)

    ## 
    pe <- estimate(list(pm,pm), lambda1 = 0,
                   data = split(d,d$gender))
    expect_equal(coef(pe),coef(e), tolerance = test.tolerance)

    ## equivariance
    pe <- estimate(list(pm,pm), lambda1 = 0, equivariance = TRUE,
                   data = split(d,d$gender))
    coef(pe)
    expect_equal(coef(pe),coef(e), tolerance = test.tolerance)
})

test_that("full penalization", {
    pm <- penalize(m, Vridge = FALSE)
    
    pe <- estimate(list(pm,pm), lambda1 = 1e8,
                   data = split(d,d$gender))
    expect_equal(names(coef0(pe)),c("1@Y~X1","1@Y~group","2@Y~X1","2@Y~group"))

    pe <- estimate(list(pm,pm), lambda1 = 1e8, equivariance = TRUE,
                   data = split(d,d$gender))
    expect_equal(names(coef0(pe)),c("1@Y~X1","1@Y~group","2@Y~X1","2@Y~group"))
})

# }}}

estimate(list(m,m), list(d,d),
         control = list(iter.max = 20, trace = 1))$opt$iterations

## 0:     1560.0775: 0.147587 0.147587  0.00000  0.00000 -0.246706  0.00000  0.00000 -0.246706
##   1:     685.16840: 0.387729 0.387729 0.0660076 0.494510 0.184659 0.0660076 0.494510 0.184659
##   2:     405.43116: 0.603980 0.603980 0.226584 0.984397 0.604636 0.226584 0.984397 0.604636
##   3:     356.57986: 0.585487 0.585487 0.667389  1.40658 0.959609 0.667389  1.40658 0.959609
##   4:     353.74740: 0.295171 0.295171 0.842373  1.60347 0.838261 0.842373  1.60347 0.838261
##   5:     353.73574: 0.295171 0.295171 0.842373  1.60347 0.822834 0.842373  1.60347 0.822834
##   6:     353.73574: 0.295171 0.295171 0.842373  1.60347 0.822952 0.842373  1.60347 0.822952
## 7:     353.73574: 0.295171 0.295171 0.842373  1.60347 0.822952 0.842373  1.60347 0.822952

estimate(list(m,m), list(d,d), control = list(iter.max = 200))$opt$iterations

#----------------------------------------------------------------------
### test-lassoRegression_proxGrad.R ends here
