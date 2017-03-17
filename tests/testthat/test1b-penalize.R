### test-Lasso_regression.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 10 2017 (14:23) 
## Version: 
## last-updated: mar 14 2017 (17:34) 
##           By: Brice Ozenne
##     Update #: 26
#----------------------------------------------------------------------
## 
### Commentary: 
## BEFORE RUNNING THE FILE:
## library(butils.base) ; package.source("lava.penalty") ;
##
## CONTENT
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
    plvm.model <- penalize(lvm.model)
    expect_equal(unname(penalty(plvm.model,no.ridge=TRUE)[["link"]]),
                 paste0("Y~X",1:12)
                 )

    plvm.model <- penalize(lvm.model, c("Y~X3","Y~X5"))
    expect_equal(unname(penalty(plvm.model,no.ridge=TRUE)[["link"]]),
                 paste0("Y~X",c(3,5))
                 )
    
    lava.options(symbols = c("Ø","*"))
    plvm.model <- penalize(lvm.model)
    expect_equal(unname(penalty(plvm.model,no.ridge=TRUE)[["link"]]),
                 paste0("YØX",1:12)
                 )

    expect_equal(penalty(plvm.model,no.ridge=TRUE,no.lasso=TRUE),
                 NULL)
    plvm.model
    lava.options(symbols = c("~","~~"))
    
    })
# }}}

# {{{ Penalize - linear regression - categorical variable

lvmCAT.model <- lvm(formula.lvm)
categorical(lvmCAT.model, K = 3, labels = letters[1:3]) <- "X1"
categorical(lvmCAT.model, K = 2, labels = letters[1:2]) <- "X2"

test_that("penalization using categorical variables",{
    plvmCAT.model <- penalize(lvmCAT.model)
    suppressWarnings(
        e <- estimate(lvmCAT.model, sim(lvmCAT.model,1e2))
    )
   
   expect_equal(unname(penalty(plvmCAT.model,no.lasso=TRUE,no.ridge=TRUE)[["link"]]),
                setdiff(names(coef(e)),coef(lvmCAT.model))
                )

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
    expect_equal(unname(unique(penalty(plvm.model)$link)),
                 c(coefReg.model, coefCov.model)
                 )
    plvm.model
    
    lava.options(symbols = c("Ø","*"))
    lvm.model <- lvm(c(Y1,Y2,Y3,Y4)~eta + X1)
    regression(lvm.model) <- Y2 ~ X2
    covariance(lvm.model) <- Y1 ~ Y2
    plvm.model <- penalize(lvm.model, covariance = TRUE, variance = TRUE, latent = TRUE)
    expect_equal(unname(unique(penalty(plvm.model)$link)),
                 c(gsub("~","Ø",coefReg.model), gsub("~~","*",coefCov.model))
                 )
    plvm.model
    lava.options(symbols = c("~","~~"))
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
