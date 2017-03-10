### test-Lasso_regression.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 10 2017 (14:23) 
## Version: 
## last-updated: mar 10 2017 (16:05) 
##           By: Brice Ozenne
##     Update #: 17
#----------------------------------------------------------------------
## 
### Commentary: 
## Test the methods:
## penalize: attribute a penalty term in a latent variable model
## penalty: extract the penalty term from a latent variable model
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(testthat)
library(lava.penalty)

# library(butils.base)
# package.source("lava.penalty", RorderDescription = FALSE)

lava.options(symbols = c("~","~~"))
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:12), collapse = "+")))
lvm.model <- lvm(formula.lvm)

# {{{ Elastic Net penalty 

# {{{ Equivalent ways to penalize a model
test_that("user interface for penalize", {
    plvm.model.1 <- penalize(lvm.model)
    plvm.model.2 <- penalize(lvm.model, value = coefReg(lvm.model, value = TRUE))
    expect_equal(plvm.model.1, plvm.model.2)
    plvm.model.3 <- penalize(lvm.model, value = Y ~ X1 + X2 + X5)
    plvm.model.4 <- penalize(lvm.model, value = list(Y ~ X1, Y ~ X2 + X5))
    expect_equal(plvm.model.3, plvm.model.4)
})
# }}}

# {{{ Extract penalty - linear regression - continous variable

test_that("no penalty", {
    expect_equal(penalty(lvm.model),
                 NULL)
    expect_equal(penalty(lvm.model, type = NULL),
                 NULL)
})

test_that("penalization using different symbols - continuous variables",{
    lava.options(symbols = c("~","~~"))
    plvm.model <- penalize(lvm.model)
    expect_equal(unname(penalty(plvm.model, type = "link")),
                 paste0("Y~X",1:12)
                 )

    plvm.model <- penalize(lvm.model, c("Y~X3","Y~X5"))
    expect_equal(unname(penalty(plvm.model, type = "link")),
                 paste0("Y~X",c(3,5))
                 )
    
    lava.options(symbols = c("Ø","*"))
    plvm.model <- penalize(lvm.model)
    expect_equal(unname(penalty(plvm.model, type = "link")),
                 paste0("YØX",1:12)
                 )

    expect_equal(penalty(plvm.model, type = "group"),
                 logical(0))
    plvm.model    
    })
# }}}

# {{{ Penalize - linear regression - categorical variable

lvmCAT.model <- lvm(formula.lvm)
categorical(lvmCAT.model, K = 3) <- "X1"
categorical(lvmCAT.model, K = 2) <- "X2"

test_that("penalization using categorical variables",{
    plvmCAT.model <- penalize(lvmCAT.model)
    
    expect_equal(unname(penalty(plvmCAT.model)),
                 c(paste0("Y",lava.options()$symbols[1],"X",3:12),
                   coefExtra(lvmCAT.model,value = TRUE))
                 )
    expect_equal(unname(penalty(plvmCAT.model, no.group = TRUE)),
                 paste0("Y",lava.options()$symbols[1],"X",3:12))
    expect_equal(unname(penalty(plvmCAT.model, no.elasticNet = TRUE)),
                 coefExtra(lvmCAT.model, value = TRUE))
    penalty(plvmCAT.model, type = c("link","group"), group = 1)
    plvmCAT.model
})
# }}}

# {{{ Penalize - LVM - continous variable 
lava.options(symbols = c("~","~~"))
lvm.model <- lvm(c(Y1,Y2,Y3,Y4)~eta + X1)
regression(lvm.model) <- Y2 ~ X2
covariance(lvm.model) <- Y1 ~ Y2


coef.model <- coef(lvm.model)
coefCov.model <- unname(coefCov(lvm.model, value = TRUE, keep.var= TRUE))
coefReg.model <- unname(coefReg(lvm.model, value = TRUE))

test_that("penalization using different symbols - continuous variables",{
    lvm.model <- lvm(c(Y1,Y2,Y3,Y4)~eta + X1)
    regression(lvm.model) <- Y2 ~ X2
    covariance(lvm.model) <- Y1 ~ Y2
    plvm.model <- penalize(lvm.model, covariance = TRUE, variance = TRUE, latent = TRUE)
    expect_equal(unname(penalty(plvm.model)),
                 c(coefReg.model, coefCov.model)
                 )
    plvm.model
    
    lava.options(symbols = c("Ø","*"))
    lvm.model <- lvm(c(Y1,Y2,Y3,Y4)~eta + X1)
    regression(lvm.model) <- Y2 ~ X2
    covariance(lvm.model) <- Y1 ~ Y2
    plvm.model <- penalize(lvm.model, covariance = TRUE, variance = TRUE, latent = TRUE)
    expect_equal(unname(penalty(plvm.model)),
                 c(gsub("~","Ø",coefReg.model), gsub("~~","*",coefCov.model))
                 )
    plvm.model
    })
# }}}

# }}}


# {{{ Nuclear norm penalty
lava.options(symbols = c("~","~~"))
coords <- expand.grid(x = 1:4, y = 1:3)

# {{{ Penalize linear regression
test_that("extract nuclear norm penalty", {
    lvm.model <- lvm(formula.lvm)
    M.coef <- matrix(coefReg(lvm.model, value = TRUE),
                     nrow = max(coords$x), ncol = max(coords$y))

    lvm.modelNuclear1 <- lvm(Y ~ X0)
    penalizeNuclear(lvm.modelNuclear1) <- M.coef
    lvm.modelNuclear1

    lvm.modelNuclear2 <- lvm(Y ~ X0)
    penalizeNuclear(lvm.modelNuclear2, coords = coords) <- coefReg(lvm.model, value = TRUE)
    lvm.modelNuclear2

    lvm.modelNuclear3 <- lvm(Y ~ X0)
    penalizeNuclear(lvm.modelNuclear3, coords = coords) <- as.formula(paste0(endogenous(lvm.model),"~",paste(exogenous(lvm.model),collapse=" + ")))
    lvm.modelNuclear3

    expect_equal(lvm.modelNuclear1$penaltyNuclear, lvm.modelNuclear2$penaltyNuclear)
    expect_equal(lvm.modelNuclear1$penaltyNuclear, lvm.modelNuclear3$penaltyNuclear)
})

# }}}

# }}}

lava.options(symbols = c("~","~~"))

#----------------------------------------------------------------------
### test-penalize.R ends here
