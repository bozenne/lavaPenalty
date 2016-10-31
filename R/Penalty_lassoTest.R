#' @title Significance test for the lasso
#'
#' @description Significance test while performing variable selection with the lasso
#' 
#' @param x a penalized latent variable model contain
#'
#' @details It is essentially a copy of the \code{covTest} function of the covTest package.
#'
#' @references 
#' R. Lockhart, J. Taylor, R. J. Tibshirani and R. Tibshirani. A significance test for the lasso. Ann Stat. 2014 42(2):413-468
#'
#' @examples 
#' m <- lvm(Y ~ X1+X2+X3+X4+X5)
#' 
#' ## simulation
#' set.seed(10)
#' mSim <- m
#' regression(mSim, Y~X1+X2+X3+X4+X5) <- list(1,0,0,0.1,0.5)
#' df <- sim(mSim,50)
#' dfs <- as.data.frame(scale(df))
#' 
#' ## lvm
#' pm <- penalize(m)
#' res <- estimate(pm, data = dfs, regularizationPath = TRUE)
#' plot(res, type = "path")
#' getPath(res)
#' lassoTest(res)
#' 
#' penalized.PathL1 <- penalized(Y ~  ., data = dfs, steps = "Park", trace = FALSE)
#' seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
#' n.lambda <- length(seq_lambda)
#' 
#' if(require(covTest)){
#'  m.lars <- lars(y = dfs[,1], x = as.matrix(dfs[,-1]), type = "lasso")
#'  a <- covTest(fitobj = m.lars, x = as.matrix(dfs[,-1]), y = dfs[,1])
#'  # covTestLasso(beta = m.lars$beta, lambda = m.lars$lambda, Y = dfs[,1], X = as.matrix(dfs[,-1]))
#'  # one difference with the LVM version is that p doesn't count the intercept
#' }
#' 
#'
#' @export
`lassoTest` <-
  function(x,...) UseMethod("lassoTest")


lassoTest.plvmfit <- function(x, ...){
  
  ## check if is a linear regression
  test.regression <- (length(latent(x)) == 0) && (length(endogenous(x)) == 1)
  if(!test.regression){
    stop("lassoTest only available for linear regression \n")
  }
  
  ## prepare argument
  index.change <- unlist(getPath(x, name = "indexChange"))
  order.lambda <- order(getPath(x)$lambda1.abs, decreasing = TRUE)
  order.lambda <- order.lambda[cumsum(!is.na(index.change[order.lambda]))>0] # cumsum(...) removes the first NA used to initialise the algorithm when increasing = FALSE
  seq.lambda <- getPath(x)$lambda1.abs[order.lambda]
  n.lambda <- length(seq.lambda)
 
  index.beta <- c(x$model$index$parBelongsTo$mean,x$model$index$parBelongsTo$reg)
  beta.all <- as.matrix(getPath(x, lambda = NULL)[order.lambda,index.beta,drop = FALSE])
  X.all <- as.matrix(x$data$model.frame[,x$model$exogenous,drop=FALSE])
  if(length(x$model$index$parBelongsTo$mean)>0){ # add intercept
    X.all <- cbind(1,X.all)
  }
  Y <- x$data$model.frame[,endogenous(x)]
  
  test <- covTestLasso(beta = beta.all, lambda = seq.lambda, Y = Y, X = X.all)
  
  ## export
  index.nameBeta <- c(NA,index.change[order.lambda])
  test <- cbind(param = colnames(beta.all)[index.nameBeta[-length(index.nameBeta)]], test)
  
  return(test)
}


covTestLasso <- function(beta, lambda, Y, X){
  
  n.lambda <- length(lambda)
  
  Tk <- data.frame(matrix(NA, nrow = n.lambda, ncol = 5))
  names(Tk) <- c("lambda","covT","cov0","statistic","p.value")
  Tk$lambda <- lambda
  
  for(iterKnot in 2:n.lambda){
    
    ## true prediction
    iter.beta <- beta[iterKnot,]
    Tk[iterKnot,"covT"] <- sum(Y  * (X %*% iter.beta))
    
    if(iterKnot > 2){ ## prediction removing following variables
      indexN0 <- which(beta[iterKnot-1,]!=0)
      
      ## linear extrapolation
      newbeta <- beta[iterKnot-1,] + (lambda[iterKnot]-lambda[iterKnot-1]) * (beta[iterKnot-1,]-beta[iterKnot-2,])/(lambda[iterKnot-1]-lambda[iterKnot-2])
      
      Tk[iterKnot,"cov0"] <- sum(Y  * (X %*% newbeta))
    }else{
      Tk[iterKnot,"cov0"] <- 0
    }
  }
  
  #### compute the statistic
  p <- NCOL(X)
  n <- NROW(Y)
  test.intercept <- apply(X,2, function(x){length(unique(x))})
  lmXY <- lsfit(x = X, y = Y, intercept = all(test.intercept!=1))
  sigma <- sqrt(sum(lmXY$res^2)/(n - p))
  
  Tk$statistic <- ((Tk$covT - Tk$cov0)/sigma^2)
  Tk$p.value <- 1 - pf(Tk$statistic, 2, n - p)
  
  return(Tk)
}