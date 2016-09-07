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
#' @param estimator the method used to compute the likelihood, gradient and hessian.
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
                          lambda1, lambda2, lambdaN, adaptive = FALSE, 
                          control = list(), estimator = "penalized", 
                          regularizationPath = FALSE, resolution_lambda1 = c(1e-1,1e-3), increasing = TRUE, reversible = FALSE, stopLambda = NULL, stopParam = NULL, exportAllPath = FALSE, 
                          method.proxGrad = "ISTA", step = 1, BT.n = 100, BT.eta = 0.8, force.descent = FALSE,
                          fit = "BIC",
                          fixSigma = FALSE, ...) {
  # names.coef <- coef(x)
  # n.coef <- length(names.coef)
  
  #### prepare control
  if("iter.max" %in% names(control) == FALSE){
    control$iter.max <- 1000
  }
  if("trace" %in% names(control) == FALSE){
    control$trace <- 0
  }else{
    control$trace <- control$trace - 1
  }
  if("constrain" %in% names(control) == FALSE){
    control$constrain <- FALSE
  }
  
  #### prepare data (scaling)
  if(control$trace>=0){cat("Scale and center dataset \n")}
  resData <- prepareData.plvm(x, data = data)
  
  x <- resData$lvm
  data <- resData$data # non scale data with factors converted to binary variables
  scaledData <- resData$scaledData # scaled data with factors converted to binary variables
  penalty <- x$penalty
  names.coef <- coef(x)
  n.coef <- length(names.coef)
  
  #### prepare penalty
  ## elastic net
  if(!missing(lambda1)){
    penalty$lambda1 <- as.numeric(lambda1)
  }
  if(!missing(lambda2)){
    penalty$lambda2 <- as.numeric(lambda2)
  }
  if(!missing(adaptive)){
    penalty$adaptive <- as.numeric(adaptive)
  }
  
  penalty$names.varCoef <- names.coef[x$index$parBelongsTo$cov]
  penalty$names.varCoef <- intersect(penalty$names.varCoef,
                                     paste(c(endogenous(x), latent(x)), 
                                           c(endogenous(x), latent(x)), 
                                           sep = ","))
  control$penalty <- penalty
  
  ## penaltyNuclear
  if(!is.null(x$penaltyNuclear$name.coef)){
    if(!missing(lambdaN)){
      x$penaltyNuclear$lambdaN <- as.numeric(lambdaN)
    }
    if(length(endogenous(x))>1){
      stop("nuclear norm only implemented for a linear model \n")
    }
    
    x$penaltyNuclear$objective <- function(coef){
      x$penaltyNuclear$FCTobjective(coef, 
                                    Y = as.vector(scaledData[,x$penaltyNuclear$name.Y,drop=TRUE]), 
                                    X = as.matrix(scaledData[,c(exogenous(x),x$penaltyNuclear$name.X),drop=FALSE]))}
    x$penaltyNuclear$gradient <- function(coef){
      x$penaltyNuclear$FCTgradient(coef, 
                                   Y = as.vector(scaledData[,x$penaltyNuclear$name.Y,drop=TRUE]), 
                                   X = as.matrix(scaledData[,c(exogenous(x),x$penaltyNuclear$name.X),drop=FALSE]))}
    control$penaltyNuclear <- x$penaltyNuclear
  }
  
  #### control/optimisation parameters 
  control$proxGrad <- list(method = method.proxGrad,
                           step = step,
                           BT.n = BT.n,
                           BT.eta = BT.eta,
                           expX = if(control$constrain){penalty$names.varCoef}else{NULL},
                           force.descent = force.descent,
                           trace = if(regularizationPath == 0){control$trace}else{FALSE},
                           fixSigma = fixSigma,
                           envir = environment() # pass x and data in case where fixSigma = TRUE
  )
  
  control$regPath <- list(type = regularizationPath,
                          resolution_lambda1 = resolution_lambda1,
                          increasing = increasing,
                          reversible = reversible,
                          stopParam = stopParam,
                          stopLambda = stopLambda,
                          fixSigma = fixSigma,
                          exportAllPath = exportAllPath,
                          trace = if(regularizationPath > 0){control$trace}else{FALSE})
  
  #### initialization
  if(is.null(control$start) || regularizationPath>0){
    if(control$trace>=0){cat("Initialization: ")}
    res.init  <- initializer.plvm(x, data = scaledData, ...)
    
    if(regularizationPath == 0){
      if(control$trace>=0){cat(" LVM where all penalized coefficient are shrinked to 0 \n")}
      control$start <- res.init$lambdaMax
    }else {
      control$regPath$beta_lambdaMax <- res.init$lambdaMax
      control$regPath$beta_lambda0 <- res.init$lambda0
      control$start <- res.init$lambdaMax
    }
  }
  
  #### define proximal operator 
  resOperator <- init.proxOperator(lambda1 = control$penalty$lambda1,  # will be NULL if control$penalty does not exist
                                   lambda2 = control$penalty$lambda2, 
                                   lambdaN = control$penaltyNuclear$lambdaN,  # will be NULL if control$penaltyNuclear does not exist
                                   group.coef = control$penalty$group.coef,
                                   regularizationPath = regularizationPath)
  control$proxOperator <- resOperator$proxOperator
  control$objectivePenalty <- resOperator$objectivePenalty
  
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
    res <- lava:::estimate.lvm(x = x, data = scaledData, 
                               method = if(regularizationPath == 0){"optim.regLL"}else{"optim.regPath"}, 
                               control = control, estimator = estimator, ...)
  }
  
  #### add elements to object
  penalty <-  control$penalty[c("name.coef", "group.coef", "lambda1", "lambda2")]
  regPath <- res$opt$regPath
  res$opt$regPath <- NULL
  
  if(regularizationPath == 0){
    if(fixSigma){
      penalty$lambda1.abs <- penalty$lambda1
      penalty$lambda2.abs <- penalty$lambda2
      penalty$lambda1 <- penalty$lambda1/sum(coef(res)[control$penalty$names.varCoef])
      penalty$lambda2 <- penalty$lambda2/sum(coef(res)[control$penalty$names.varCoef])
    }else{
      penalty$lambda1.abs <- penalty$lambda1*sum(coef(res)[control$penalty$names.varCoef])
      penalty$lambda2.abs <- penalty$lambda2*sum(coef(res)[control$penalty$names.varCoef])
    }
  }
  
  #### estimate the best model according to the fit parameter
  if(regularizationPath > 0 && !is.null(fit)){
    if(control$trace>=0){cat("Best penalized model according to the",fit,"criteria",if(control$trace>=1){"\n"})}
    resLambda <- calcLambda(path = regPath, model = x, fit = fit, data.fit = data, trace = control$trace)
    res <- resLambda$optimum$lvm
    resLambda$optimum$lvm <- NULL
    regPath <- resLambda
    if(control$trace>=0){cat(" - done \n")}
  }
  
  #### export
  class(res) <- append("plvmfit", class(res))
  res$regularizationPath <- regPath
  res$penalty <- penalty
  res$penaltyNuclear <- control$penaltyNuclear
  res$x <- x
  
  return(res)
}


#' @title Intialization of the parameter of the plvm
#'
#' @description Compute the coefficients of the non and the completely regularized latent variable model
#' 
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u+x1+x2+x3
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1 ~ y2
#' 
#' pm <- penalize(m)
#' initializer(pm, data = sim(pm, 1e2))
#' 
initializer.plvm <- function(x, data, ...){
  
  names.coef <- coef(x)
  n.coef <- length(names.coef)
  n.data <- NROW(data)
  penalty <- x$penalty
  
  #### normal model
  if(n.data > n.coef){
    suppressWarnings(
      initLVM <- try(lava:::estimate.lvm(x = x, data = data, ...), silent = TRUE)
    )
    
    if(("try-error" %in% class(initLVM) == FALSE)){ # should also check convergence
      start_lambda0 <- setNames(rep(0,n.coef),names.coef)
      start_lambda0[names(coef(initLVM))] <- coef(initLVM)
    }else{ 
      start_lambda0 <- NULL
    }
  }else{
    start_lambda0 <- NULL
  }
  
  #### hight dimensional model
  x0 <- x
  start_lambdaMax <- setNames(rep(0,n.coef),names.coef)
  
  ## removed penalized variables
  for(iter_link in penalty$name.coef){
    x0 <- rmLink(x0, iter_link)
  }
  
  ## estimate the model
  suppressWarnings(
    x0.fit <- lava:::estimate.lvm(x = x0, data = data, ...)
  )
  newCoef <- coef(x0.fit, level = 9)[,"Estimate"]
  newCoef <- newCoef[names(newCoef) %in% names(start_lambdaMax)] # only keep relevant parameters
  start_lambdaMax[names(newCoef)] <- newCoef
  # start_lambdaMax[] <- coef(x0.fit, level = 9)[,"Estimate"]
  
  return(list(lambda0 = start_lambda0,
              lambdaMax = start_lambdaMax))
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
  
  if(any(manifest(x) %in% names(data) == FALSE)){
    stop("prepareData.lvm: arguments \'data\' and \'x\' are incompatible \n",
         "variables: ",paste(manifest(x)[manifest(x) %in% names(data) == FALSE], collapse = " ")," not found in \'data\' \n")
  }
  
  #### convert categorical variables to dummy variables
  resC2D <- lava:::categorical2dummy(x, data)
  
  index.numeric <- intersect(manifest(x), manifest(resC2D$x))
  indexOld.factor <- setdiff(manifest(x),  manifest(resC2D$x))
  indexNew.factor <- setdiff(manifest(resC2D$x), manifest(x))
  
  test.factor <- length(indexNew.factor)>0
  
  if(test.factor){
    if(any(endogenous(x) %in% indexOld.factor == TRUE)){
      stop("prepareData.lvm: endogenous variables must not be categorical \n",
           "incorrect variables: ",paste(endogenous(x)[endogenous(x) %in% indexOld.factor == TRUE], collapse = " "),"\n")
    }
    ls.factor <- lapply(indexOld.factor, function(var){unique(data[,var])})
    names(ls.factor) <- indexOld.factor
    conversion.factor <- sapply(indexNew.factor, renameFactor, ls.factor)
  }else{
    conversion.factor <- NULL
  }
  
  x <- resC2D$x
  data <- resC2D$data
  
  #### rescale data
  if(class(data)[1] != "data.frame"){data <- as.data.frame(data)}
  
  if(length(index.numeric)>0){
    value.center <- sapply(index.numeric, function(x){do.call(method.center,args = list(data[[x]]))})
    value.scale <- sapply(index.numeric, function(x){do.call(method.scale,args = list(data[[x]]))})
    data[, index.numeric] <- scale(data[, index.numeric, drop = FALSE], center = value.center, scale = value.scale)
  }
  
  #### update the penalty according to the dummy variables
  if(test.factor){
    penalty <- x$penalty
    
    ## find the new coefficients to penalize
    name.Newlinks <- coef(x)
    
    OldLinksPenalty.factors <- setdiff(penalty$name.coef, name.Newlinks)
    NewExogePenalty.factor <- lapply(1:length(OldLinksPenalty.factors), function(x){
      formulaTempo <- as.formula(OldLinksPenalty.factors[x])
      newVar <- names(conversion.factor)[conversion.factor %in% all.vars(formulaTempo)[2]]
      return(paste0(all.vars(formulaTempo)[1],"~",newVar))
    })  
    
    name.NewlinksPenalty <- c(intersect(penalty$name.coef, name.Newlinks),
                              unlist(NewExogePenalty.factor))
    
    x$penalty$V <- NULL
    res <- penalize(x, value = name.NewlinksPenalty)
    
    x$penalty <- res$penalty
  }
  
  #### export
  return(list(data = resC2D$data,
              scaledData = data, 
              scale = value.scale,
              center = value.center,
              lvm = x))
}



orthoData_glmPath <- function(model, name.Y, allCoef, penaltyCoef, data){
  
  ## function
  extractVar <- function(names){
    names.formula <- grep("~", names ,fixed = TRUE)  
    if(length(names.formula)>0){
      names[names.formula] <- unlist(lapply(strsplit(names[names.formula], split = "~", fixed = TRUE),"[",2))
    }
    
    return(names)
  }
  
  ## preparation
  n <- nrow(data)
  n.coef <- length(allCoef)
  names.interceptCoef <- intersect(allCoef,coef(model)[model$index$parBelongsTo$mean])
  n.interceptCoef <- length(names.interceptCoef)
  names.covCoef <- intersect(allCoef,coef(model)[model$index$parBelongsTo$cov])
  if(length(model$latent)>0){
    names.latentCoef <- allCoef[sapply(names(model$latent), grep, x = allCoef, fixed = TRUE)]
  }else{
    names.latentCoef <- NULL
  }
  
  var.penalized <-  extractVar(penaltyCoef)
  var.unpenalized <- setdiff(extractVar(setdiff(allCoef,c(names.covCoef,names.latentCoef))), 
                             c(var.penalized))
  
  ## rebuild data
  X_tempo <- data[,setdiff(c(var.penalized,var.unpenalized),names.interceptCoef), drop = FALSE]
  
  if(n.interceptCoef>0){
    X_tempo <- cbind(matrix(1, nrow = n, ncol = n.interceptCoef),
                     X_tempo)
  }
  names(X_tempo)[1:n.interceptCoef] <- names.interceptCoef
  
  ## distinguish penalized from non penalized
  penalized <-  as.matrix(X_tempo[,var.penalized, drop = FALSE])
  unpenalized <-  as.matrix(X_tempo[,var.unpenalized, drop = FALSE])
  
  ## orthogonlize
  orthogonalizer <- solve(crossprod(unpenalized), crossprod(unpenalized, penalized))
  penalized <- penalized - unpenalized %*% orthogonalizer
  
  ## starting coefficients
  mu.X <- setNames(rep(0,n.coef), allCoef)
  lm.fitted <- lm.fit(y = data[[name.Y]], x = unpenalized)
  mu.X[names(mu.X) %in% penaltyCoef == FALSE] <- c(coef(lm.fitted), var(lm.fitted$residuals)) ### issue with the latent variable here!!
  
  ## scale
  sd.X <- setNames(rep(1,n.coef), allCoef)
  
  varNI.penalized <- setdiff(var.penalized, names.interceptCoef)
  index.penalized <- setdiff(which(names(mu.X) %in% penaltyCoef),
                             which(names(mu.X) %in% names.interceptCoef) )
  
  if(length(varNI.penalized)>0){
    sd.X[index.penalized] <- sqrt(apply(penalized[,varNI.penalized, drop = FALSE], 2, var)*(n-1)/n)
    penalized <- sweep(penalized[,varNI.penalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.penalized])
  }
  
  varNI.unpenalized <- setdiff(var.unpenalized, names.interceptCoef)
  index.unpenalized <- setdiff(which(names(mu.X) %in% penaltyCoef == FALSE), 
                               which(names(mu.X) %in% c(names.interceptCoef, names.covCoef))
  )
  
  if(length(varNI.unpenalized)>0){
    sd.X[index.unpenalized] <- sqrt(apply(unpenalized[,varNI.unpenalized, drop = FALSE], 2, var)*(n-1)/n)
    unpenalized <- sweep(unpenalized[,varNI.unpenalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.unpenalized])
  }
  mu.X <- mu.X * sd.X
  
  ## update initial dataset
  data[,var.penalized] <- penalized
  if(length(setdiff(var.unpenalized, names.interceptCoef))>0){
    data[,setdiff(var.unpenalized, names.interceptCoef)] <- unpenalized
  }
  
  ## lambda
  lambda1 <- setNames(rep(0,n.coef), allCoef)
  lambda1[which(names(mu.X) %in% penaltyCoef)] <- 1/sd.X[which(names(mu.X) %in% penaltyCoef)]
  lambda2 <- setNames(rep(0,n.coef), allCoef)
  lambda2[which(names(mu.X) %in% penaltyCoef)] <- 1/(sd.X[which(names(mu.X) %in% penaltyCoef)]^2)
  
  ## export
  return(list(data = data,
              orthogonalizer = orthogonalizer,
              sd.X = sd.X,
              mu.X = mu.X,
              lambda1 = lambda1,
              lambda2 = lambda2))
}
