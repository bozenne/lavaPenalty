### Penalty_hook.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 10 2017 (17:00) 
## Version: 
## last-updated: mar 10 2017 (15:21) 
##           By: Brice Ozenne
##     Update #: 169
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


# {{{ lava.penalty.estimate.hook
#' @title Hook to estimate a penalized lvm model
#' @description Add LP to dataset, update the estimator for handling LP, and find initialisation.
#' 
#' @examples 
#' 
#' m <- lvm()
#' m <- regression(m,y='y1',x='x'%++%1:2)
#' m <- regression(m,y='y1',x='z'%++%1:5)
#' 
#' # simul
#' set.seed(10)
#' d <- sim(m,5e2)
#' d <- as.data.frame(scale(d))
#' 
#' suppressWarnings(
#' start <- coef(estimate(m, data = d, control = list(iter.max = 0)))
#' )
#' 
#' # penalized lvm 
#' mp <- penalize(m)
#' system.time(
#' lassoLVM <- estimate(mp, data = d, lambda1 = 10, control = list(trace = 2, start = start[coef(mp)], iter.max = 1), estimator = "gaussian2")
#' )
#' 
#' 
#' # reduced penalized lvm
#' mp.red <- reduce(mp)
#' system.time(
#' lassoRLVM <- estimate(mp.red, data = d, lambda1 = 10, control = list(trace = 2, start = start[coef(mp.red)], iter.max = 1))
#' )
#' 
#' coef(lassoRLVM)-coef(lassoLVM)[names(coef(lassoRLVM))]
#' lassoRLVM$opt$iterations
#' 
#' @export
lava.penalty.estimate.hook <- function(x,data,weight,weight2,estimator,...) {

    dots <- list(...)
    if("plvm" %in% class(x) && ("proxGrad" %in% names(dots$optim) || "regPath" %in% names(dots$optim)) ){ # 
        test.regularizationPath <- !is.null(dots$optim$regPath)
        start <- dots$optim$start
        trace <- dots$optim$trace

        #### initialization
        if(is.null(start) || test.regularizationPath){
            if(trace>0){cat("Initialization: ")}
            res.init <- initializer.plvm(x, data = data, regularizationPath = test.regularizationPath,
                                         constrain.variance = dots$optim$proxGrad$constrain.variance,
                                         name.variance = dots$optim$proxGrad$name.variance,
                                         ...)
            if(test.regularizationPath){
                dots$optim$regPath$start.lambda0 <- res.init$lambda0
                dots$optim$regPath$start.lambdaMax <- res.init$lambdaMax
            }else {
                if(trace>0){cat(" LVM where all penalized coefficient are shrinked to 0 \n")}
                dots$optim$start <- res.init$lambdaMax
            }
        }
    
        #### add penalty in control
        dots$optim$penalty <- x$penalty
        dots$optim$penaltyNuclear <- x$penaltyNuclear
    }

    return(c(list(x=x,
                  data=data,
                  weight=weight,
                  weight2=weight2,
                  estimator=estimator),dots))   
}
# }}}

# {{{ lava.penalty.post.hook

#' @export
lava.penalty.post.hook <- function(x){
    
    if(any(c("proximal gradient","EPSODE") %in% x$opt)){

        penalty <- x$control$penalty
        penaltyNuclear <- x$control$penaltyNuclear

        trace <- x$control$trace
        constrain.lambda <- x$control$proxGrad$constrain.lambda
        name.variance <- x$control$proxGrad$name.variance
        test.regularizationPath <- !is.null(x$control$regPath)
        fit <- x$control$regPath$fit
        regPath <- x$opt$regPath

    
    #### add elements to object
    
    if(test.regularizationPath){
            
      if(!is.null(fit)){ #### estimate the best model according to the fit parameter
        
        if(trace>=0){cat("Best penalized model according to the",fit,"criteria",if(x$control$trace>=1){"\n"})}
        resLambda <- calcLambda(path = regPath,
                                fit = fit,
                                model = x$model,
                                data.fit = x$data$model.frame,
                                trace = trace+1)
        res <- resLambda$optimum$lvm
        resLambda$optimum$lvm <- NULL
        regPath <- resLambda
        if(trace>=0){cat(" - done \n")}
        
      }
      
    }else{
      if(constrain.lambda){
        penalty$lambda1.abs <- penalty$lambda1
        penalty$lambda2.abs <- penalty$lambda2
        penalty$lambda1 <- penalty$lambda1/sum(coef(x)[name.variance])
        penalty$lambda2 <- penalty$lambda2/sum(coef(x)[name.variance])
      }else{
        penalty$lambda1.abs <- penalty$lambda1*sum(coef(x)[name.variance])
        penalty$lambda2.abs <- penalty$lambda2*sum(coef(x)[name.variance])
      }
    }
    
        #### export
        x$regularizationPath <- regPath
        x$penalty <- penalty
        x$penaltyNuclear <- penaltyNuclear
        class(x) <- append("plvmfit", class(x))   
    }
  
    return(x)
  
}

# }}}



#----------------------------------------------------------------------
### Penalty_hook.R ends here
