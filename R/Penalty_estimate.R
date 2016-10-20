`initializer` <-
  function(x,...) UseMethod("initializer")

#' @title Penalized lvm model
#' @aliases estimate estimate.plvm
#' @description Estimate a penalized lvm model
#'
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param lambda1 lasso penalization parameter
#' @param lambda2 ridge penalization parameter
#' @param lambdaN Nuclear norm penalization parameter
#' @param adaptive should the coefficient of the adaptive lasso be returned instead of the coefficient of the lasso
#' @param control control/optimization parameters
#' @param regularizationPath should the regularization path be computed. If so the argument lambda1 is ignored but not lambda2.
#' @param resolution_lambda1,increasing,stopLambda,stopParam arguments for the EPSODE function (see \code{EPSODE})
#' @param method.proxGrad,step,BT.n,BT.eta,force.descent arguments for the proximal gradient algorithm (see \code{proxGrad}) 
#' @param fit criterion to decide of the optimal model to retain among the penalized models.
#' @param fixSigma should the variance parameter be fixed at 1 ? Only works for regression models. [temporary]
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
#' examples/EX_Penalty_estimate.R
#' 
#' @export
estimate.plvm <- function(x, data, 
                          lambda1 = NULL, lambda2 = NULL, lambdaN = NULL, adaptive = FALSE, 
                          control = list(), 
                          regularizationPath = FALSE, resolution_lambda1 = lava.options()$EPSODE$resolution_lambda1, increasing = TRUE,stopLambda = NULL, stopParam = NULL, 
                          method.proxGrad = lava.options()$proxGrad$method,
                          fit = lava.options()$calcLambda$fit,
                          fixSigma = FALSE, ...) {
  
  #### prepare penalty
  # update penalty (for new variables)
  newV <- penalty(x, type = "V")
  penalty(x, type = "V") <- newV[coef(x),coef(x)]
  
  # update with potential user input
  penalty <- penalty(x, type = c("lambda1","lambda2","adaptive"))
  if(!is.null(lambda1)){penalty$lambda1 <- as.numeric(lambda1)}
  if(!is.null(lambda2)){penalty$lambda2 <- as.numeric(lambda2)}
  if(adaptive){penalty$adaptive <- as.numeric(adaptive)}
  
  ## penaltyNuclear
  penaltyNuclear <- penalty(x, type = "lambdaN", nuclear = TRUE, keep.list = TRUE)
  if(!is.null(lambdaN)){penaltyNuclear$lambdaN <- as.numeric(lambdaN)}
  
  #### prepare control 
  # overwrite lava standard setting for the optimisation
  controlUser <- control
  control <- list(iter.max = lava.options()$proxGrad$iter.max,
                  constrain = lava.options()$constrain,
                  trace = lava.options()$trace,
                  start = setNames(rep(NA, length(coef(x))), coef(x)))
  if (length(controlUser) > 0) {
    control[names(controlUser)] <- controlUser
  }
  control$trace <- control$trace - 1 
  
  control$penalty <- penalty
  control$penaltyNuclear <- penaltyNuclear
  
  control$proxGrad <- list(expX = if(control$constrain){penalty$names.varCoef}else{NULL},
                           trace = if(regularizationPath == 0){control$trace}else{FALSE},
                           method = method.proxGrad,
                           fixSigma = fixSigma,
                           envir = environment() # pass x and data in case where fixSigma = TRUE
  )
  
  if(regularizationPath){
    
    # if linear regression then use LARS algorithm
    test.regression <- (length(endogenous(x)) == 1) && (length(latent(x)) == 0)
    if(test.regression && "resolution_lambda1" %in% names(match.call()) == FALSE){
      resolution_lambda1 <- c(1,1e-10)
      fixSigma <- TRUE
      increasing <- FALSE
    }
    
    control$regPath <- list(resolution_lambda1 = resolution_lambda1,
                            increasing = increasing,
                            stopParam = stopParam,
                            stopLambda = stopLambda,
                            fixSigma = fixSigma,
                            fit = fit,
                            trace = control$trace)
  }
  
  #### prepare data (scaling, convert character to dummy variables)
  if(control$trace>=0){cat("Scale and center dataset \n")}
  resData <- prepareData.plvm(x, data = data)
  x <- resData$lvm
  data <- resData$data
  
  #### optimisation
  if(all(c("objective", "gradient", "hessian") %in% names(list(...)))){
    
    if(regularizationPath == 0){
      res <- optim.regLL(start = control$start, 
                         objective = list(...)$objective, gradient = list(...)$gradient, hessian = list(...)$hessian, 
                         control = control)
    }else{
      res <- optim.regPath(objective = list(...)$objective, gradient = list(...)$gradient, hessian = list(...)$hessian, 
                           control = control)
    }
    return(res)
    
  }else{
    res <- estimate.lvm(x = x, data = data, 
                        method = if(regularizationPath == 0){"optim.regLL"}else{"optim.regPath"}, 
                        control = control, ...)
  }
  
  #### export
  return(res)
}




#' @title Prepare the data for the estimate function
#'  
#' @description Scale the data, update the penalty term according to the presence of factors. 
#' 
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param method.center function used to center the data
#' @param method.scale function used to scale the data
#' 
prepareData.plvm <- function(x, data, method.center = "mean", method.scale = "sd"){
  
  if(any(manifest(x, lp = FALSE) %in% names(data) == FALSE)){
    stop("prepareData.lvm: arguments \'data\' and \'x\' are incompatible \n",
         "variables: ",paste(manifest(x)[manifest(x) %in% names(data) == FALSE], collapse = " ")," not found in \'data\' \n")
  }
  
  #### convert categorical variables to dummy variables
  resC2D <- lava:::categorical2dummy(x, data)
  
  index.numeric <- intersect(manifest(x, lp = FALSE), manifest(resC2D$x, lp = FALSE))
  indexOld.factor <- setdiff(manifest(x, lp = FALSE),  manifest(resC2D$x, lp = FALSE))
  indexNew.factor <- setdiff(manifest(resC2D$x, lp = FALSE), manifest(x, lp = FALSE))
  
  test.factor <- length(indexNew.factor)>0
  
  if(test.factor){
    if(any(endogenous(x) %in% indexOld.factor == TRUE)){
      stop("prepareData.lvm: endogenous variables must not be categorical \n",
           "incorrect variables: ",paste(endogenous(x)[endogenous(x) %in% indexOld.factor == TRUE], collapse = " "),"\n")
    }
    ls.factor <- lapply(indexOld.factor, function(var){unique(data[[var]])})
    names(ls.factor) <- indexOld.factor
    conversion.factor <- sapply(indexNew.factor, renameFactor, ls.level = ls.factor)
  }else{
    conversion.factor <- NULL
  }
  
  x <- resC2D$x
  data <- resC2D$data
  
  #### rescale data
  if(class(data)[1] != "data.frame"){data <- as.data.frame(data)}
  if(length(index.numeric)>0){
    value.center <- sapply(index.numeric, function(x){do.call(method.center,args = list(na.omit(data[[x]])))})
    value.scale <- sapply(index.numeric, function(x){do.call(method.scale,args = list(na.omit(data[[x]])))})
    data[, index.numeric] <- scale(data[, index.numeric, drop = FALSE], center = value.center, scale = value.scale)
  }
  
  #### update the penalty according to the dummy variables
  if(test.factor){
    link.penalty <- penalty(x, type = "link")
    group.penalty <- penalty(x, type = "group")
    
    ## find the new coefficients to penalize
    name.Newlinks <- coef(x)
    
    OldLinksPenalty.factors <- setdiff(link.penalty, name.Newlinks)
    NewExogePenalty.factor <- lapply(1:length(OldLinksPenalty.factors), function(x){
      formulaTempo <- as.formula(OldLinksPenalty.factors[x])
      newVar <- names(conversion.factor)[conversion.factor %in% all.vars(formulaTempo)[2]]
      return(paste0(all.vars(formulaTempo)[1],"~",newVar))
    })  
    names(NewExogePenalty.factor) <- OldLinksPenalty.factors
    
    ## redefine the group of penalty
    name.commonPenalty <- intersect(link.penalty, name.Newlinks)
    name.NewlinksPenalty <- c(name.commonPenalty,
                              unname(unlist(NewExogePenalty.factor)))
    group.NewlinksPenalty <- setNames(rep(NA, length(name.NewlinksPenalty)),
                                      name.NewlinksPenalty)
    group.NewlinksPenalty[name.commonPenalty] <- group.penalty[na.omit(match(link.penalty,name.commonPenalty))]
    
    for(iterG in OldLinksPenalty.factors){
      groupG <- group.penalty[match(iterG,link.penalty)]
      if(groupG>1){ # if already a group, keep it
        group.NewlinksPenalty[NewExogePenalty.factor[[iterG]]] <- groupG
      }else{ # else create a new group
        group.NewlinksPenalty[NewExogePenalty.factor[[iterG]]] <- floor(max(c(na.omit(group.NewlinksPenalty), group.penalty))) + 1
      }
    }
    
    ## redefine V matrix
    penalty(x, type = "V") <- NULL
    res <- penalize(x, value = name.NewlinksPenalty, group = group.NewlinksPenalty, add = FALSE)
    
    penalty(x, type = NULL) <- penalty(res, type = NULL)
  }
  
  #### export
  return(list(data = data,
              scale = value.scale,
              center = value.center,
              lvm = x))
}


