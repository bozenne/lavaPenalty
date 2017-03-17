### test-getVars.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar 13 2017 (10:15) 
## Version: 
## last-updated: mar 14 2017 (12:59) 
##           By: Brice Ozenne
##     Update #: 17
#----------------------------------------------------------------------
## 
### Commentary:
## BEFORE RUNNING THE FILE:
## library(butils.base) ; package.source("lava.penalty") ;
## path <- file.path(butils.base::path_gitHub(),"lava.penalty","tests")
## source(file.path(path,"FCT.R"))
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(testthat)
library(data.table)

context("#### test extractors #### \n")

m <- lvm(Y ~ X1 + X2 + X3 + X4)

test_that("all continuous", {
    res <- getIvar.lvm(m)
    test.type <- res[coefReg(m),type]
    expect_true(all(test.type=="continuous"))
})

m1 <- m
test_that("one categorical - numeric labels", {
    categorical(m1, K = 3, labels = 1:3) <- ~X1
    res <- getIvar.lvm(m1)
    expect_true(all(res[link %in% setdiff(coefReg(m1, value = TRUE),"Y~X1"),type]=="continuous"))
    expect_true(all(res[externalLink %in% coefExtra(m1, value = TRUE),type]=="categorical"))
    expect_false("Y~X1" %in% res[["link"]])

    suppressWarnings(em1 <- estimate(m1, sim(m1,1e2)))
    a <- res[["link"]]
    b <- setdiff(names(coef(em1)),coefExtra(em1,value=TRUE))
    expect_equal(sort(a),sort(b))
    
    res <- getIvar.lvm(m1, link = "Y~X1")
    expect_false("Y~X1" %in% res[["link"]])
    expect_true(all(res[link %in% coefExtra(m1, value = TRUE),externalLink]=="categorical"))
})

m2 <- m
test_that("two categorical - numeric labels", {
    categorical(m2, K = 3, labels = 2:4) <- ~X1
    categorical(m2, K = 2, labels = letters[1:2]) <- ~X3
    res <- getIvar.lvm(m2)
    expect_false("Y~X1" %in% res[["link"]])
    expect_false("Y~X3" %in% res[["link"]])
    expect_true(all(res[link %in% coefExtra(m2, value = TRUE),externalLink]=="categorical"))

    suppressWarnings(em2 <- estimate(m2, sim(m2,1e2)))
    a <- res[["link"]]
    b <- setdiff(names(coef(em2)),coefExtra(em2,value=TRUE))
    expect_equal(sort(a),sort(b))

    m2.bis <- lava_categorical2dummy(m2, sim(m2, 1))$x
    res <- getIvar.lvm(m2.bis, data = sim(m2, 1))
})

#----------------------------------------------------------------------
### test-getVars.R ends here
