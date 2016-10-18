#' @title Hook to estimate a penalized lvm model
#' @description Add LP to dataset, update the estimator for handling LP, and find initialisation.
#' 
#' @examples 
#' 
#' m <- lvm()
#' m <- regression(m,y='y1',x='x'%++%1:2)
#' m <- regression(m,y='y1',x='z'%++%1:50)
#' 
#' # simul
#' set.seed(10)
#' d <- sim(m,5e2)
#' d <- as.data.frame(scale(d))
#' 
#' # penalized lvm 
#' mp <- penalize(m)
#' system.time(
#' lassoLVM <- estimate(mp, data = d, lambda1 = 10, control = list(trace = 1, iter.max = 1000, constrain = TRUE))
#' )
#' 
#' # reduced penalized lvm
#' mp.red <- reduce(mp)
#' system.time(
#' lassoRLVM <- estimate(mp.red, data = d, lambda1 = 10, control = list(trace = 2, iter.max = 1000, constrain = TRUE))
#' )
#' @export
lava.penalty.estimate.hook <- function(x,data,weight,weight2,estimator,...) {
  
  dots <- list(...)
  
  if("plvm" %in% class(x) && "penalty" %in% names(dots$optim) ){
    
    control <- dots$optim
    penalty <- control$penalty
    penaltyNuclear <- control$penaltyNuclear
    proxGrad <- control$proxGrad
    regPath <- control$regPath
    
    test.regularizationPath <- !is.null(regPath)
    
    # update the object
    penalty(x, type = "lambda1") <- penalty$lambda1
    penalty(x, type = "lambda2") <- penalty$lambda2
    penalty(x, type = "adaptive") <- penalty$adaptive
    penalty(x, nuclear = TRUE, type = "lambdaN") <- penaltyNuclear$lambdaN
    control[["penalty"]] <- NULL
    
    names.coef <- coef(x)
    n.coef <- length(names.coef)
    
    #### initialization
    if(all(is.na(control$start)) || test.regularizationPath){
      if(control$trace>=0){cat("Initialization: ")}
      res.init  <- initializer.plvm(x, data = data, ...)
      
      if(test.regularizationPath){
        control$regPath$beta_lambdaMax <- res.init$lambdaMax
        control$regPath$beta_lambda0 <- res.init$lambda0
        control$start <- res.init$lambdaMax
      }else {
        if(control$trace>=0){cat(" LVM where all penalized coefficient are shrinked to 0 \n")}
        control$start <- res.init$lambdaMax
      }
    }
    
    #### define proximal operator 
    resOperator <- init.proxOperator(lambda1 = control$penalty$lambda1,  # will be NULL if control$penalty does not exist
                                     lambda2 = control$penalty$lambda2, 
                                     lambdaN = control$penaltyNuclear$lambdaN,  # will be NULL if control$penaltyNuclear does not exist
                                     group.coef = control$penalty$group.coef,
                                     regularizationPath = test.regularizationPath)
    
    penalty(x, type = "proxOperator") <- resOperator$proxOperator
    penalty(x, type = "objectivePenalty") <- resOperator$objectivePenalty
    
    #### update
    dots$optim <- control
    dots$optim$penalty <- penalty(x, type = NULL)
    dots$optim$penaltyNuclear <- penalty(x, nuclear = TRUE, type = NULL)
    dots$optim$proxGrad <- proxGrad
    dots$optim$regPath <- regPath
  }
  
  return(c(list(x=x,data=data,weight=weight,weight2=weight2,estimator=estimator),dots)) 
  
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
  
  ## intialization 
  suppressWarnings(
    start_init <- lava:::estimate.lvm(x = x, data = data, control = list(iter.max = 0))
  )
  
  names.coef  <- names(coef(start_init))
  n.coef <- length(names.coef)
  n.data <- NROW(data)
  
  #### normal model
  
  if(n.data > n.coef){
    suppressWarnings(
      initLVM <- try(lava:::estimate.lvm(x = x, data = data, ...), silent = TRUE)
    )
    
    if(("try-error" %in% class(initLVM) == FALSE)){ # should also check convergence
      start_lambda0 <- coef(initLVM)
    }else{ 
      start_lambda0 <- NULL
    }
  }else{ 
    start_lambda0 <- NULL
  }
  
  #### high dimensional model
  x0 <- x
  start_lambdaMax <- setNames(rep(0,n.coef),names.coef)
  
  ## removed penalized variables
  for(iter_link in penalty(x, type = "link")){
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


lava.penalty.post.hook <- function(x){
  
  if("penalty" %in% names(x$control)){
    penalty <- x$control$penalty
    penaltyNuclear <- x$control$penaltyNuclear
    proxGrad <- x$control$proxGrad
    regPath <- x$control$regPath
    
    test.regularizationPath <- !is.null(regPath)
    
    #### add elements to object
    
    if(test.regularizationPath){
      regPath <- res$opt$regPath
      
      if(!is.null(regPath$fit)){ #### estimate the best model according to the fit parameter
        
        if(x$control$trace>=0){cat("Best penalized model according to the",fit,"criteria",if(x$control$trace>=1){"\n"})}
        resLambda <- calcLambda(path = regPath, fit = regPath$fit, model = x, data.fit = eval(x$call$data), trace = x$control$trace+1)
        res <- resLambda$optimum$lvm
        resLambda$optimum$lvm <- NULL
        regPath <- resLambda
        if(x$control$trace>=0){cat(" - done \n")}
        
      }
      
    }else{
      if(proxGrad$fixSigma){
        penalty$lambda1.abs <- penalty$lambda1
        penalty$lambda2.abs <- penalty$lambda2
        penalty$lambda1 <- penalty$lambda1/sum(coef(x)[penalty$var.coef])
        penalty$lambda2 <- penalty$lambda2/sum(coef(x)[penalty$var.coef])
      }else{
        penalty$lambda1.abs <- penalty$lambda1*sum(coef(x)[penalty$var.coef])
        penalty$lambda2.abs <- penalty$lambda2*sum(coef(x)[penalty$var.coef])
      }
      
    }
    
    #### export
    class(x) <- append("plvmfit", class(x))
    x$regularizationPath <- regPath
    x$penalty <- penalty
    x$penaltyNuclear <- x$control$penaltyNuclear
  }
  
  return(x)
  
}

