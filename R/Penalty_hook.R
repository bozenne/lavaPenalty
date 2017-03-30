### Penalty_hook.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 10 2017 (17:00) 
## Version: 
## last-updated: mar 30 2017 (16:42) 
##           By: Brice Ozenne
##     Update #: 223
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
lava.penalty.estimate.hook <- function(x,data,...) {
    dots <- list(...)
    if("plvm" %in% class(x) && ("proxGrad" %in% names(dots$optim) || "regPath" %in% names(dots$optim)) ){ #
        control.regPath <- dots$optim$regPath
        test.regularizationPath <- !is.null(control.regPath) && (control.regPath$algorithm == "EPSODE") 
        increasing <- dots$optim$regPath$increasing
        start <- dots$optim$start
        trace <- dots$optim$trace

        #### initialization
        if(is.null(start) || test.regularizationPath){
            if(trace>0){cat("Initialization: ")}
            dots$optim$start <- initialize.start(x, data = data,
                                                 regularizationPath = test.regularizationPath, increasing = increasing,
                                                 constrain.variance = dots$optim$proxGrad$constrain.variance,
                                                 name.variance = dots$optim$proxGrad$name.variance,
                                                 ...)
            if(trace>0){cat("done \n")}
        }

        #### add penalty in control
        dots$optim$penalty <- x$penalty
        dots$optim$penaltyNuclear <- x$penaltyNuclear
    }

    return(c(list(x=x, data=data),dots))   
}
# }}}

# {{{ lava.penalty.post.hook

#' @export
lava.penalty.post.hook <- function(x){

    if(any(c("EPSODE","proximal gradient") %in% x$opt$algorithm )){

        #### update x
        x$regularizationPath <- x$opt$regPath
        x$penalty <- x$control$penalty
        x$penaltyNuclear <- x$control$penaltyNuclear
        class(x) <- append("plvmfit", class(x))   

        x$opt$regPath <- NULL
        x$control$penalty <- NULL
        x$control$penaltyNuclear <- NULL

        #### add elements to object
        trace <- x$control$trace
        constrain.lambda <- x$control$proxGrad$constrain.lambda
        name.variance <- x$control$proxGrad$name.variance
        fit <- x$control$regPath$fit

        if(is.path(x)){
            ## estimate the best model according to the fit parameter
            if(!is.null(fit)){
                x <- calcLambda(x, fit = fit, trace = trace)                                
            }
      
        }else{
            if(constrain.lambda){
                x$penalty$lambda1.abs <- x$penalty$lambda1
                x$penalty$lambda2.abs <- x$penalty$lambda2
                x$penalty$lambda1 <- x$penalty$lambda1/sum(coef(x)[name.variance])
                x$penalty$lambda2 <- x$penalty$lambda2/sum(coef(x)[name.variance])
            }else{
                x$penalty$lambda1.abs <- x$penalty$lambda1*sum(coef(x)[name.variance])
                x$penalty$lambda2.abs <- x$penalty$lambda2*sum(coef(x)[name.variance])
            }
        }
    
        #### export
    }
  
    return(x)
  
}

# }}}



#----------------------------------------------------------------------
### Penalty_hook.R ends here
