
# {{{ optim.regLL
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
    
    ## only keep the relevant elements in control
    control.proxGrad <- control$proxGrad
    penalty <- control$penalty
    penaltyNuclear <- control$penaltyNuclear
    
    # {{{ update penalty according to the coefficients present in start
    name.coef <- names(start)
    n.coef <- length(name.coef)
    newPenalty <- updatePenaltyVariable.penaltyL12(penalty, name.coef = name.coef)

    index.tempo <- setNames(penalty(x, type = "group"), name.groupLasso)
        index.penaltyG <- tapply(index.tempo,index.tempo, function(x){
            list(setNames(match(names(x), name.coef), names(x)))
        })
    
    newPenaltyNuclear <- initializer.penaltyNuclear(penaltyNuclear, name.coef = name.coef)        
    if(!identical(newPenaltyNuclear$lambdaN,0)){
        hessian <- NULL
    }
    browser()
    # }}}
    
    # {{{ Proximal gradient 
    ## generate proximal operator
    resInit <- initializeOperator(lambda1 = newPenalty$lambda1,
                                  lambda2 = newPenalty$lambda2,
                                  lambdaG = newPenalty$lambdaG, index.penaltyG = newPenalty$index.penaltyG,
                                  lambdaN = newPenaltyNuclear$lambdaN, index.penaltyN = newPenaltyNuclear$index.penaltyN, 
                                  nrow = newPenaltyNuclear$nrow, ncol = newPenaltyNuclear$ncol,
                                  constrain.lambda = control$proxGrad$constrain.lambda,
                                  constrain.variance = control$proxGrad$constrain.variance,
                                  index.variance = control$proxGrad$name.variance)
    # resInit$objectivePenalty(start)
    # resInit$proxOperator(start, step = 0.1)
    
    ## run descent
    res <- proxGrad(start = start, proxOperator = resInit$proxOperator,
                    hessian = hessian, gradient = gradient, objective = objective, 
                    control = control$proxGrad)
    # }}}

    # {{{ Adaptive proximal gradient
    if(penalty$adaptive){
        ## generate proximal operator
        resInit <- initializeOperator(lambda1 = newPenalty$lambda1/abs(res$par),
                                      lambda2 = newPenalty$lambda2,
                                      lambdaG = newPenalty$lambdaG, index.penaltyG = newPenalty$index.penaltyG,
                                      lambdaN = newPenaltyNuclear$lambdaN, index.penaltyN = newPenaltyNuclear$index.penaltyN, 
                                      nrow = newPenaltyNuclear$nrow, ncol = newPenaltyNuclear$ncol,
                                      constrain.lambda = control$proxGrad$fixSigma,
                                      constrain.variance = control$proxGrad$constrain,
                                      index.variance = control$proxGrad$name.variance)

        ## run descent    
        res <- proxGrad(start = res$par, proxOperator = resInit$proxOperator,
                        hessian = hessian, gradient = gradient, objective = objective, 
                        control = control$proxGrad)
    }

    # }}}

    ## update objective with penalty
    res$objectivePenalty <- resInit$objectivePenalty
    
    ## export
    #attr(res$message,"par") <- res$par
    #res$par <- res$par[name.coef]
    return(res)
}
# }}}

# {{{ optim.regPath

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
  penalty <- control$penalty
  penaltyNuclear <- control$penaltyNuclear
  regPath <- control$regPath
  
  #### update the penalty according to start 
  # (some coefficient have been removed as they are chosen as a reference)
  newPenalty <- initPenalty(name.coef = name.coef, pen = penalty, regPath = regPath, penNuclear = penaltyNuclear)
    if(!is.null(regPath)){ ## to be modified
        regPath$beta_lambda0 <- regPath$beta_lambda0[name.coef, drop = FALSE]
        regPath$beta_lambdaMax <- regPath$beta_lambdaMax[name.coef, drop = FALSE]
    }
  
  ##
  penalty <- newPenalty$penalty
  regPath <- newPenalty$regPath
  PGcontrols <- c("iter.max","trace", "constrain", "proxGrad", "proxOperator", "regPath")
  control <- control[names(control) %in% PGcontrols]
  
  ## nuisance parameter
  if(regPath$fixSigma){
    #if(length(penalty$names.varCoef)>1){stop("Cannot fixSigma when having several variance parameters \n")}
    indexNuisance <- which(name.coef %in% penalty$var.coef)
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
    newPenalty <- initPenalty(start = start, pen = penalty, penNuclear = penaltyNuclear)
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
    
    
    
  }
  
  ## EPSODE
  link.penalty <- penalty(penalty, type = "link")
  indexPenalty <-  which(name.coef %in% link.penalty)
  
  resEPSODE <- EPSODE(beta_lambda0 = regPath$beta_lambda0, beta_lambdaMax = regPath$beta_lambdaMax,
                      objective = objective, gradient = gradient, hessian = hessian, 
                      V = penalty$V, indexPenalty = indexPenalty, indexNuisance = indexNuisance, lambda2 = penalty$lambda2,
                      resolution_lambda1 = regPath$resolution_lambda1, increasing = regPath$increasing, stopLambda = regPath$stopLambda, stopParam = regPath$stopParam,
                      constrain = control$constrain, trace = control$trace)
  
  ## estimation of the nuisance parameter and update lambda/lambda.abs
  if(length(indexNuisance)>0){
    if(control$trace>=0){cat("Estimation of the nuisance parameter ")}
    resEPSODE$path[,name.coef[indexNuisance]] <- sapply(1:NROW(resEPSODE$path), function(x){
      optim.Nuisance.plvm(x = control$proxGrad$envir$x, data = control$proxGrad$envir$data,
                          coefEstimated = unlist(resEPSODE$path[x,name.coef]),
                          coefNuisance = name.coef[indexNuisance], coefPenalty = penalty(penalty, type = "link"),
                          control = control, ...)}
    )
    if(control$trace>=0){cat("- done \n")}
    resEPSODE$path$lambda1 <- resEPSODE$path$lambda1.abs/resEPSODE$path[,name.coef[indexNuisance]]
    resEPSODE$path$lambda2 <- resEPSODE$path$lambda2.abs/resEPSODE$path[,name.coef[indexNuisance]]
  }else{
    sumSigma <- apply(resEPSODE$path[,penalty$var.coef,drop=FALSE],1,sum)
    resEPSODE$path$lambda1.abs <- resEPSODE$path$lambda1*sumSigma
    resEPSODE$path$lambda2.abs <- resEPSODE$path$lambda2*sumSigma
  }
  
  #### conversion to regPath object
  message <- resEPSODE$message
  regPath <- list(path = resEPSODE$path,
                  increasing = regPath$increasing, 
                  lambda = if(regPath$fixSigma){"lambda1.abs"}else{"lambda1"},
                  penCoef = link.penalty,
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

# }}}

# {{{ optim.Nuisance.plvm
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
optim.Nuisance.plvm <- function(x, data, 
                                coefEstimated, coefNuisance, coefPenalty,
                                control, ...){
  
  control$trace <- FALSE
  names.coef <- setdiff(names(coefEstimated), coefNuisance)
  n.coef <- length(names.coef)
  
  fit.coef <- coefEstimated[names.coef, drop = FALSE]
  
  
  xConstrain <- x
  class(xConstrain) <- setdiff(class(x),"plvm")
  
  index.keep <- which(setdiff(names.coef, coefPenalty) %in% names.coef)[1]  # leave the first coefficient to the model
  if(length(index.keep) == 0){
    warning("optim.Nuisance: all coefficients are penalized - may lead to incorrect estimation of the nuisance parameter \n")
  }
  
    for(iter_p in setdiff(1:n.coef,index.keep)){
        cat("[TODO] optim.Nuisance.plvm to update \n")
        browser()
    if(names(fit.coef)[iter_p] %in% coefPenalty && fit.coef[iter_p]==0){
      xConstrain <- rmLink(xConstrain, var1 = names.coef[iter_p], simplify = TRUE)
    }else{
      ## remove the link from the linear predictor (if any)
      if("lvm.reduced" %in% class(xConstrain) && names.coef[iter_p] %in% lp(xConstrain, type = "link")){
        xConstrain <- cancelLP(xConstrain, link = names.coef[iter_p])
        xConstrain <- addLink(xConstrain, var1 = names.coef[iter_p], covariance = FALSE)
      }
      xConstrain <- setLink(xConstrain, var1 = names.coef[iter_p], value = fit.coef[iter_p])
      
    }
  }
  
  elvm <- estimate(xConstrain, data = data, control = control)
  
  return(coef(elvm)[coefNuisance])
  
}
# }}}



