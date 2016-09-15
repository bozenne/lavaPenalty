#' @title Proximal gradient algorithm
#' 
#' @description Estimate a penalized lvm model using a proximal gradient algorithm 
#' 
#' @param start -exported by estimate.lvm
#' @param objective -exported by estimate.lvm
#' @param gradient -exported by estimate.lvm
#' @param hessian -exported by estimate.lvm
#' @param control -exported by estimate.lvm 
#' @param ... -exported by estimate.lvm
#' 
#' @export 
optim.regLL <- function(start, objective, gradient, hessian, control, ...){
  PGcontrols <- c("iter.max","trace", "constrain", "proxOperator", "objectivePenalty",
                  "proxGrad")
  
  penalty <- control$penalty
  penaltyNuclear <- control$penaltyNuclear
  control <- control[names(control) %in% PGcontrols]
  
  name.coef <- names(start)
  n.coef <- length(name.coef)
  
 #### Proximal gradient algorithm
  # update penalty
  newPenalty <- initPenalty(start = start, penalty = penalty, penaltyNuclear = penaltyNuclear)
  
  # force equivariant
  if(control$proxGrad$fixSigma){
    index.constrain <- which(name.coef %in% penalty$names.varCoef)
    if(control$trace>=0){cat("constrains: lambda1=lambda1/sum(", paste(name.coef[index.constrain], collapse = " "),")\n")}
  }else{
    index.constrain <- NULL
  }
  
  
  proxOperator <- function(x, step){
    control$proxOperator(x, step = step,
                         lambdaN = newPenalty$lambdaN, lambda1 = newPenalty$lambda1, lambda2 = newPenalty$lambda2, 
                         test.penaltyN = newPenalty$test.penaltyN, test.penalty1 = newPenalty$test.penalty1, test.penalty2 = newPenalty$test.penalty2,
                         nrow = penaltyNuclear$nrow, ncol = penaltyNuclear$ncol,
                         index.constrain = index.constrain, type.constrain = control$constrain, expX = control$proxGrad$expX)
  }
  if(!is.na(newPenalty$lambdaN)){
    hessian <- penaltyNuclear$hessian
    gradient <- penaltyNuclear$gradient
    objective <- penaltyNuclear$objective
    if(length(index.constrain)>0){index.constrain <- which(names(newPenalty$start) %in% name.coef[index.constrain])}
  }
  
  if(control$trace>=0){cat("Proximal gradient ")}
  
    res <- proxGrad(start = newPenalty$start, proxOperator = proxOperator, method = control$proxGrad$method,
                  hessian = hessian, gradient = gradient, objective = objective, 
                  iter.max = control$iter.max, trace = control$proxGrad$trace)
  if(control$trace>=0){cat("- done \n")}

  if(penalty$adaptive){
    
    proxOperator <- function(x, step){
      control$proxOperator(x, step = step,
                           lambdaN = newPenalty$lambdaN, lambda1 = newPenalty$lambda1/abs(res$par), lambda2 = newPenalty$lambda2, 
                           test.penaltyN = newPenalty$test.penaltyN, test.penalty1 = newPenalty$test.penalty1, test.penalty2 = newPenalty$test.penalty2,
                           index.constrain = index.constrain, type.constrain = control$constrain, expX = control$proxGrad$expX)
    }
    
    if(control$trace>=0){cat("Proximal gradient (adaptive) ")}
    res <- proxGrad(start = res$par, proxOperator = proxOperator, method = control$proxGrad$method,
                    hessian = hessian, gradient = gradient, objective = objective,
                    iter.max = control$iter.max, trace = control$proxGrad$trace)
    if(control$trace>=0){cat("- done \n")}
  }

  ## update objective with penalty - one value for each penalty
  objective.pen <- lapply(control$objectivePenalty, 
                          function(fct){do.call(fct, args = list(res$par, 
                                                                 lambdaN = newPenalty$lambdaN, lambda1 = newPenalty$lambda1, lambda2 = newPenalty$lambda2, 
                                                                 test.penalty1 = newPenalty$test.penalty1, test.penalty2 = newPenalty$test.penalty2, test.penaltyN = newPenalty$test.penaltyN,
                                                                 nrow = penaltyNuclear$nrow, ncol = penaltyNuclear$ncol,
                                                                 expX = control$proxGrad$expX)
                          )})
  
  res$objective <- objective(res$par) + sum(unlist(objective.pen))

  ### export
  attr(res$message,"par") <- res$par
  res$par <- res$par[name.coef]
  
  return(res)
}

#' @title Regularization path
#' 
#' @description Estimate the regularization path associated to a LVM
#' 
#' @param start -exported by estimate.lvm
#' @param objective -exported by estimate.lvm
#' @param gradient -exported by estimate.lvm
#' @param hessian -exported by estimate.lvm
#' @param control -exported by estimate.lvm
#' @param ... -exported by estimate.lvm
#' 
#' @export
optim.regPath <- function(start, objective, gradient, hessian, control, ...){
  
  name.coef <- names(control$regPath$beta_lambdaMax)
  penaltyNuclear <- control$penaltyNuclear
  
  #### update the penalty according to start 
  # (some coefficient have been removed as they are chosen as a reference)
  newPenalty <- initPenalty(name.coef = name.coef, penalty = control$penalty, regPath = control$regPath, penaltyNuclear = NULL)
 
  ##
  penalty <- newPenalty$penalty
  regPath <- newPenalty$regPath
  PGcontrols <- c("iter.max","trace", "constrain", "proxGrad", "proxOperator", "regPath")
  control <- control[names(control) %in% PGcontrols]
 
  ## nuisance parameter
  if(regPath$fixSigma){
    #if(length(penalty$names.varCoef)>1){stop("Cannot fixSigma when having several variance parameters \n")}
    indexNuisance <- which(name.coef %in% penalty$names.varCoef)
  }else{
    indexNuisance <- NULL
  }
  
  ## reestimate beta
  if(penalty$lambda2 > 0){
    
    # force equivariant
    if(control$proxGrad$fixSigma){
      index.constrain <- which(name.coef %in% penalty$names.varCoef)
    }else{
      index.constrain <- NULL
    }
    
    penalty$lambda1 <- 0
    newPenalty <- initPenalty(start = start, penalty = penalty, penaltyNuclear = NULL)
    proxOperator <- function(x, step){
      control$proxOperator(x, step = step,
                           lambdaN = newPenalty$lambdaN, lambda1 = newPenalty$lambda1, lambda2 = newPenalty$lambda2, 
                           test.penaltyN = newPenalty$test.penaltyN, test.penalty1 = newPenalty$test.penalty1, test.penalty2 = newPenalty$test.penalty2,
                           nrow = penaltyNuclear$nrow, ncol = penaltyNuclear$ncol,
                           index.constrain = index.constrain, type.constrain = control$constrain, expX = control$proxGrad$expX)
    }
    
    regPath$beta_lambda0 <- do.call("proxGrad",
                                    list(start = start, proxOperator = proxOperator, method = control$proxGrad$method,
                                         hessian = hessian, gradient = gradient, objective = objective,
                                         iter.max = control$iter.max, trace = FALSE))$par
    
    penalty$lambda1 <- 1e10
    newPenalty <- initPenalty(start = start, penalty = penalty, penaltyNuclear = NULL)
    proxOperator <- function(x, step){
      control$proxOperator(x, step = step,
                           lambdaN = newPenalty$lambdaN, lambda1 = newPenalty$lambda1, lambda2 = newPenalty$lambda2, 
                           test.penaltyN = newPenalty$test.penaltyN, test.penalty1 = newPenalty$test.penalty1, test.penalty2 = newPenalty$test.penalty2,
                           nrow = penaltyNuclear$nrow, ncol = penaltyNuclear$ncol,
                           index.constrain = index.constrain, type.constrain = control$constrain, expX = control$proxGrad$expX)
    }
    regPath$beta_lambdaMax <- do.call("proxGrad",
                                      list(start = start, proxOperator = proxOperator,  method = control$proxGrad$method,
                                           hessian = hessian, gradient = gradient, objective = objective,
                                           iter.max = control$iter.max, trace = FALSE))$par
    
  }
  
  ## EPSODE
  resEPSODE <- EPSODE(beta_lambda0 = regPath$beta_lambda0, beta_lambdaMax = regPath$beta_lambdaMax,
                      objective = objective, gradient = gradient, hessian = hessian, 
                      V = penalty$V, indexPenalty = which(name.coef %in% penalty$name.coef), indexNuisance = indexNuisance, lambda2 = penalty$lambda2,
                      resolution_lambda1 = regPath$resolution_lambda1, increasing = regPath$increasing, stopLambda = regPath$stopLambda, stopParam = regPath$stopParam,
                      constrain = control$constrain, trace = control$trace)

 ## estimation of the nuisance parameter and update lambda/lambda.abs
 if(length(indexNuisance)>0){
   if(control$trace>=0){cat("Estimation of the nuisance parameter ")}
   resEPSODE$path[,name.coef[indexNuisance]] <- sapply(1:NROW(resEPSODE$path), function(x){
     optim.Nuisance.lvm(x = control$proxGrad$envir$x, data = control$proxGrad$envir$scaledData,
                        coefEstimated = unlist(resEPSODE$path[x,name.coef]),
                        coefNuisance = name.coef[indexNuisance], coefPenalty = penalty$name.coef,
                        control = control)}
   )
   if(control$trace>=0){cat("- done \n")}
   resEPSODE$path$lambda1 <- resEPSODE$path$lambda1.abs/resEPSODE$path[,name.coef[indexNuisance]]
   resEPSODE$path$lambda2 <- resEPSODE$path$lambda2.abs/resEPSODE$path[,name.coef[indexNuisance]]
 }else{
   sumSigma <- apply(resEPSODE$path[,penalty$names.varCoef,drop=FALSE],1,sum)
   resEPSODE$path$lambda1.abs <- resEPSODE$path$lambda1*sumSigma
   resEPSODE$path$lambda2.abs <- resEPSODE$path$lambda2*sumSigma
 }
  
  #### conversion to regPath object
  message <- resEPSODE$message
  regPath <- list(path = resEPSODE$path,
                  increasing = regPath$increasing, 
                  lambda = if(regPath$fixSigma){"lambda1.abs"}else{"lambda1"},
                  penCoef = penalty$name.coef,
                  performance = NULL,
                  optimum = NULL)
  class(regPath) <- "regPath"
  
  res <- list(par = control$regPath$beta_lambdaMax,
              convergence = 0,
              iterations = 0,
              evaluations = c("function" = 0, "gradient" = 0),
              message = message, 
              regPath = regPath,
              objective = NA)
  
  ### export
  return(res)
}

#' @title Estimate nuisance parameter
#' 
#' @description Use a constrained LVM to estimate nuisance parameters
#' 
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param coefEstimated the estimated mean coefficients
#' @param coefNuisance the name of the nuisance parameters to be estimated
#' @param coefPenalty the name of the penalized parameters
#' @param control control arguments to be passed to estimate.lvm
#' 
optim.Nuisance.lvm <- function(x, data, 
                              coefEstimated, coefNuisance, coefPenalty,
                              control){
  
  control$trace <- FALSE
  names.coef <- setdiff(names(coefEstimated), coefNuisance)
  n.coef <- length(names.coef)
  
  fit.coef <- coefEstimated[names.coef, drop = FALSE]
 
  
  xConstrain <- x
  class(xConstrain) <- "lvm"
  
  index.keep <- which(setdiff(names.coef, coefPenalty) == names.coef)[1]  # leave the first coefficient to the model
  if(length(index.keep) == 0){
    warning("optim.Nuisance: all coefficients are penalized - may lead to incorrect estimation of the nuisance parameter \n")
  }
  
  for(iter_p in setdiff(1:n.coef,index.keep)){
    if(names(fit.coef)[iter_p] %in% coefPenalty && fit.coef[iter_p]==0){
      xConstrain <- rmLink(xConstrain, var1 = names.coef[iter_p])
    }else{
      xConstrain <- setLink(xConstrain, var1 = names.coef[iter_p], value = fit.coef[iter_p])
    }
  }
 
  elvm <- estimate(xConstrain, data = data, control = control)
  
  return(coef(elvm)[coefNuisance])
  
}

#' @title Initialise the penalty 
#' 
#' @description Update the penalty according to the initialization of estimate
#'  
#' @param start initialization value of the parameters
#' @param name.coef names of the parameters
#' @param regPath control parameters for the regularization path
#' @param penalty definition of the penalty
#' @param penaltyNuclear definition of the nuclear norm penalty
#' 
initPenalty <- function(penalty, penaltyNuclear, start = NULL, name.coef = NULL, regPath = NULL){
  
  if(is.null(name.coef)){name.coef <- names(start)}
  n.coef <- length(name.coef)
  
  ## lasso
  if(!is.null(penalty$name.coef)){
    
    
    ## check
    if(any(penalty$name.coef %in% name.coef == FALSE)){
      warning("initPenalty: some penalty will not be applied because the corresponding parameter is used as a reference \n",
              "non-applied penalty: ",paste(setdiff(penalty$name.coef, name.coef), collapse = " "),"\n")
      penalty$group.coef <- penalty$group.coef[penalty$name.coef %in% name.coef]
      penalty$name.coef <- penalty$name.coef[penalty$name.coef %in% name.coef]
      penalty$var.coef <- penalty$var.coef[penalty$var.coef %in% name.coef] # useless
    }
    
    if(!is.null(regPath)){
      penalty$V <- penalty$V[name.coef, name.coef,drop = FALSE]
      regPath$beta_lambda0 <- regPath$beta_lambda0[name.coef, drop = FALSE]
      regPath$beta_lambdaMax <- regPath$beta_lambdaMax[name.coef, drop = FALSE]
    }
    
    lambdaN <- NA
    lambda1 <- rep(0, n.coef)
    lambda1[name.coef %in% penalty$name.coef] <- penalty$lambda1
    lambda2 <- rep(0, n.coef)
    lambda2[name.coef %in% penalty$name.coef] <- penalty$lambda2
    index.constrain <- rep(NA, n.coef) # useless
    index.constrain[name.coef %in% penalty$var.coef] <- penalty$var.coef # useless
    
    ## grouped lasso: set lasso indexes to 0
    test.penaltyN <- NULL
    test.penalty1 <- setNames(rep(0, n.coef), name.coef)
    test.penalty1[penalty$name.coef] <- penalty$group.coef
    test.penalty2 <- lambda2>0
  }
  
  if(!is.null(penaltyNuclear$name.coef)){
    
    if(any(penaltyNuclear$name.coef %in% name.coef == FALSE)){
    start <- c(start[setdiff(name.coef,penalty$names.varCoef)],
               setNames(rep(0,length(penaltyNuclear$name.coef)),penaltyNuclear$name.coef),
               start[penalty$names.varCoef])  
    }
    
    name.coef <- names(start)
    n.coef <- length(start)
    
    lambdaN <- penaltyNuclear$lambdaN
    lambda1 <- NA
    lambda2 <- NA
    
    test.penaltyN <- which(name.coef %in% penaltyNuclear$name.coef)
    test.penalty1 <- NULL
    test.penalty2 <- NULL
    index.constrain <- NULL
  }
  
  return(list(start = start, regPath = regPath, penalty = penalty,
              lambdaN = lambdaN, lambda1 = lambda1, lambda2 = lambda2, index.constrain = index.constrain,
              test.penaltyN = test.penaltyN, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2))
}




