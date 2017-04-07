
# {{{ doc
#' @title Penalized lvm model
#' @name estimate
#' @aliases estimate estimate.plvm
#' @description Estimate a penalized lvm model
#'
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param estimator the method used to compute the likelihood and its first two derivatives.
#' @param lambda1 lasso penalization parameter
#' @param lambda2 ridge penalization parameter
#' @param lambdaG group lasso penalization parameter
#' @param lambdaN nuclear norm penalization parameter
#' @param adaptive should the coefficient of the adaptive lasso be returned instead of the coefficient of the lasso
#' @param regularizationPath should the regularization path be computed. If so the argument lambda1 is ignored but not lambda2.
#' @param fit the criterion used to retain the optimal lambda1 value when the regularization path is used.
#' @param control control/optimization parameters
#' @param control.proxGrad control values for the proximal gradient algorithm (see \code{proxGrad}) 
#' @param control.EPSODE control values for the EPSODE algorithm (see \code{EPSODE}) 
#' @param equivariance should the regularization parameters be multiplied by the first variance parameter and the first variance parameter fixed to 0.
#' @param ... additional arguments to be passed to lava:::estimate.lvm
#' 
#' @details Available estimators are:
#' \itemize{
#'  \item{"penalized"}{ the hessian is computed as E[t(S) S]}
#'  \item{"numDeriveSimple"}{ the hessian is computed using numerical derivatives}
#'  \item{"numDeriveRichardson"}{ the hessian is computed using numerical derivatives (Richardson method)}
#'  \item{"explicit"}{ the hessian is computed using an explicit formula}
#' }   
#' 
#' @references 
#' Zhou 2014 - A generic Path Algorithm for Regularized Statistical Estimation \cr
#' Park 2007 - L1-regularization path algorithm for GLM
#' 
#' @return a plvmfit object
#' 
#' @example 
#' R/examples/EX_Penalty_estimate.R
# }}}

# {{{ estimate.plvm
#' @rdname estimate
#' @export
estimate.plvm <- function(x, data, estimator = "gaussian1",
                          lambda1, lambda2, lambdaG, lambdaN, adaptive = FALSE, 
                          regularizationPath = FALSE, fit = lava.options()$calcLambda$fit,
                          control = list(), control.proxGrad = list(), control.EPSODE = list(), 
                          equivariance = FALSE, ...) {
  
  if(identical(regularizationPath,TRUE)){
    method <- "optim.EPSODE"
  }else if(identical(regularizationPath,FALSE)){
    method <- "optim.proxGrad"
  }else{
    method <- "optim.proxGradPath"
  }
  
  # {{{ prepare control
  controlUser <- control
  
  controlUser.proxGrad <- control.proxGrad
  nameUser <- names(controlUser.proxGrad)
  nameLava <- names(lava::lava.options()$proxGrad)
  if(any(nameUser %in% nameLava == FALSE)){
    stop("Some elements in the arguments \'control.proxGrad\' does not match the names of lava.options()$proxGrad \n",
         "elements : \"",paste(nameUser[nameUser %in% nameLava == FALSE],collapse = "\" \""),"\" \n")
  }
  
  controlUser.EPSODE <- control.EPSODE
  nameUser <- names(controlUser.EPSODE)
  nameLava <- names(lava::lava.options()$EPSODE)
  if(any(nameUser %in% nameLava == FALSE)){
    stop("Some elements in the arguments \'control.EPSODE\' does not match the names of lava.options()$EPSODE \n",
         "elements : \"",paste(nameUser[nameUser %in% nameLava == FALSE],collapse = "\" \""),"\" \n")
  }
  
  ## lava for initialisation
  if("trace" %in% names(controlUser) == FALSE){
    control$trace <- FALSE
  }
  if("constrain" %in% names(controlUser) == FALSE){
    control$constrain <- lava::lava.options()$constrain
  }
  
  ## proximal gradient
  control$proxGrad <- c(lava::lava.options()$proxGrad,
                        adaptive = adaptive,
                        equivariance = equivariance,
                        constrain.variance = control$constrain,
                        name.variance = list(coefVar(x, value = TRUE))
  )
  if (length(controlUser.proxGrad) > 0) { # update with user input
    control$proxGrad[names(controlUser.proxGrad)] <- controlUser.proxGrad
  }
  
  ## regularization path
  if(method == "optim.EPSODE"){   
    control$regPath <- c(algorithm = "EPSODE",
                         lava::lava.options()$EPSODE,
                         equivariance = equivariance,
                         constrain.variance = control$constrain,
                         name.variance = list(coefVar(x, value = TRUE)),
                         envir = environment(),
                         fit = fit)
    
    # if linear regression then use LARS algorithm
    test.regression <- (length(endogenous(x)) == 1) && (length(latent(x)) == 0)
    if(test.regression && "lars" %in% names(controlUser.EPSODE) == FALSE){
      control$regPath$lars <- TRUE
      control$regPath$equivariance <- TRUE
      control$regPath$increasing <- FALSE
    }
    if (length(control.EPSODE) > 0) { # update with user input
      control$regPath[names(control.EPSODE)] <- control.EPSODE
    }
  }else if(method == "optim.proxGradPath"){       
    control$regPath <- c(algorithm = "proxGrad",
                         grid = list(regularizationPath),
                         warmUp = lava.options()$proxGradPath$warmUp)
  }else{
    control$regPath <- NULL
  }
  # }}}
  
  
  # {{{ prepare penalty
  ### update penalty (for factor variables)
  # the original lvm and data are needed
  x <- initializeFactor.penaltyL12(x, data = data,
                                   trace = control$trace)
  table.penalty <- penalty(x, nuclear = FALSE)
  
  ### update the penalty with potential user input
  if("lasso" %in% table.penalty$penalty && method == "optim.proxGrad"){
    penalty(x, type = "lambda1", add = FALSE) <- lambda1
  }
  if("ridge" %in% table.penalty$penalty){
    penalty(x, type = "lambda2", add = FALSE) <- lambda2
  }
  if("group lasso" %in% table.penalty$penalty){
    penalty(x, type = "lambdaG", add = FALSE) <- lambdaG
  }
  if(!is.null(penalty(x, nuclear = TRUE, type = "link"))){    
    penalty(x, type = "lambdaN", nuclear = TRUE) <- lambdaN
  }
  # }}}
  
  
  # {{{ prepare data
  ## update data substituing binary indicators to factors (+ scaling)
  if(control$trace>0){cat("Scale and center dataset \n")}
  resData <- initializeData(x, data = data)
  x <- resData$lvm
  data <- resData$data    
  
  if(any(table.penalty$link %in% coef(x) == FALSE)){
    stop("unknown penalized links \n",
         "links: \"",paste(table.penalty$link[table.penalty$link %in% coef(x) == FALSE], collapse = "\" \""),"\"\n")
  }    
  # }}}
  
  
  
    # {{{ optimisation
  res <- lava_estimate.lvm(x = x, data = data, estimator = estimator,
                           method = method,
                           control = control, quick = FALSE, index = TRUE, ...)
  # }}}
  
  ## export
  return(res)    
}

# }}}

# {{{ estimate.pmultigroup
#' @rdname estimate
#' @export
estimate.pmultigroup <- function(x, data2, estimator = "gaussian1",
                                 lambda1, lambda2, lambdaG, lambdaN, adaptive = FALSE, 
                                 regularizationPath = FALSE, fit = lava.options()$calcLambda$fit,
                                 control = list(), control.proxGrad = list(), control.EPSODE = list(), 
                                 equivariance = FALSE, ...) {

  if(identical(regularizationPath,TRUE)){
    method <- "optim.EPSODE"
  }else if(identical(regularizationPath,FALSE)){
    method <- "optim.proxGrad"
  }else{
    method <- "optim.proxGradPath"
  }
  
  # {{{ prepare control
  controlUser <- control
  
  controlUser.proxGrad <- control.proxGrad
  nameUser <- names(controlUser.proxGrad)
  nameLava <- names(lava::lava.options()$proxGrad)
  if(any(nameUser %in% nameLava == FALSE)){
    stop("Some elements in the arguments \'control.proxGrad\' does not match the names of lava.options()$proxGrad \n",
         "elements : \"",paste(nameUser[nameUser %in% nameLava == FALSE],collapse = "\" \""),"\" \n")
  }
  
  controlUser.EPSODE <- control.EPSODE
  nameUser <- names(controlUser.EPSODE)
  nameLava <- names(lava::lava.options()$EPSODE)
  if(any(nameUser %in% nameLava == FALSE)){
    stop("Some elements in the arguments \'control.EPSODE\' does not match the names of lava.options()$EPSODE \n",
         "elements : \"",paste(nameUser[nameUser %in% nameLava == FALSE],collapse = "\" \""),"\" \n")
  }
  
  ## lava for initialisation
  if("trace" %in% names(controlUser) == FALSE){
    control$trace <- FALSE
  }
  if("constrain" %in% names(controlUser) == FALSE){
    control$constrain <- lava::lava.options()$constrain
  }
  
  ## proximal gradient
  control$proxGrad <- c(lava::lava.options()$proxGrad,
                        adaptive = adaptive,
                        equivariance = equivariance,
                        constrain.variance = control$constrain,
                        name.variance = list(coefVar(x, value = TRUE))
  )
  if (length(controlUser.proxGrad) > 0) { # update with user input
    control$proxGrad[names(controlUser.proxGrad)] <- controlUser.proxGrad
  }
  
  ## regularization path
  if(method == "optim.EPSODE"){   
    control$regPath <- c(algorithm = "EPSODE",
                         lava::lava.options()$EPSODE,
                         equivariance = equivariance,
                         constrain.variance = control$constrain,
                         name.variance = list(coefVar(x, value = TRUE)),
                         envir = environment(),
                         fit = fit)
    
    # if linear regression then use LARS algorithm
    test.regression <- (length(endogenous(x)) == 1) && (length(latent(x)) == 0)
    if(test.regression && "lars" %in% names(controlUser.EPSODE) == FALSE){
        control$regPath$lars <- TRUE
        control$regPath$equivariance <- TRUE
        control$regPath$increasing <- FALSE
    }
    if (length(control.EPSODE) > 0) { # update with user input
        control$regPath[names(control.EPSODE)] <- control.EPSODE
    }
  }else if(method == "optim.proxGradPath"){       
      control$regPath <- c(algorithm = "proxGrad",
                           grid = list(regularizationPath),
                           warmUp = lava.options()$proxGradPath$warmUp)
  }else{
      control$regPath <- NULL
  }
    # }}}
    
    # {{{ prepare penalty
    # note: initializeFactor is already done in the hook: lava.penalty.multigroup.hook
    table.penalty <- penalty(x, nuclear = FALSE)
  
    ### update the penalty with potential user input
    if("lasso" %in% table.penalty$penalty && method == "optim.proxGrad"){
        penalty(x, type = "lambda1", add = FALSE) <- lambda1
    }
    if("ridge" %in% table.penalty$penalty){
        penalty(x, type = "lambda2", add = FALSE) <- lambda2
    }
    if("group lasso" %in% table.penalty$penalty){
        penalty(x, type = "lambdaG", add = FALSE) <- lambdaG
    }
    if(!is.null(penalty(x, nuclear = TRUE, type = "link"))){    
        penalty(x, type = "lambdaN", nuclear = TRUE) <- lambdaN
    }

    control$penalty <- x$penalty
    control$penaltyNuclear <- initPenaltNuclear()
    # }}}
  
    
    # {{{ prepare data
    ## scale the dataset
    if(control$trace>0){cat("Scale and center dataset \n")}
    
    resData <- initializeData(x)
    x <- resData$lvm

    if(any(table.penalty$link %in% unlist(x$name) == FALSE)){
        stop("unknown penalized links \n",
             "links: \"",paste(table.penalty$link[table.penalty$link %in% unlist(x$name) == FALSE], collapse = "\" \""),"\"\n")
    }    
    # }}}
  
  
  
    # {{{ optimisation
    for(iModel in 1:length(x$lvm)){
        class(x$lvm[[iModel]]) <- setdiff(class(x$lvm[[iModel]]), "plvm")
    }
    res <- lava_estimate.multigroup(x = x, data = data2, estimator = estimator,
                                    method = method,
                                    control = control, quick = FALSE, index = TRUE, ...)
    # }}}
  
  ## export
  return(res)    
}

# }}}





