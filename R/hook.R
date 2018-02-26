### Penalty_hook.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 10 2017 (17:00) 
## Version: 
## last-updated: apr 20 2017 (13:48) 
##           By: Brice Ozenne
##     Update #: 281
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


# {{{ lavaPenalty.estimate.hook
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
lavaPenalty.estimate.hook <- function(x,data,...) {
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

# {{{ lavaPenalty.multigroup.hook

#' @title Hook to estimate a multigroup penalized lvm model
#' @description Merge the penalty of the different models into one
#' 
#' @export
lavaPenalty.multigroup.hook <- function(mg, x, ...) {
    dots <- list(...)
    
    ls.penaltyNuclear <- lapply(x, penalty, nuclear = TRUE)
    if(any(unlist(lapply(ls.penaltyNuclear, is.null)==FALSE))){    
      stop("nuclear penalty for multigroup LVM not implemented")
    }else{
      mg$penaltyNuclear <- list(lambdaN = 0,
                                nrow = NULL,
                                ncol = NULL,
                                name.reduce = NULL,
                                endogenous = NULL,
                                link = NULL)
    }
    
    ls.penalty <- lapply(x, penalty)
    if(any(unlist(lapply(ls.penalty, is.null)==FALSE))){
        penalty.mg <- list(lambda1 = numeric(0),
                           lambda2 = numeric(0),
                           lambdaG = numeric(0),
                           Vlasso = NULL,
                           Vridge = NULL,
                           Vgroup = NULL)
        class(penalty.mg) <- "penaltyL12"

        ### update penalty (for factor variables)
        n.model <- length(x)
        for(iModel in 1:n.model){ # iModel <- 1
            mg$lvm[[iModel]]$penalty <- initializeFactor.penaltyL12(x[[iModel]],
                                                                    data = mg$data[[iModel]],
                                                                    trace = FALSE)
        }
        ls.penalty <- lapply(mg$lvm, penalty)

        ##
        if( "lasso" %in% unlist(lapply(ls.penalty,"[[","penalty")) ){
            penalty.mg$Vlasso <- initVcoef.multigroup(mg, type = "lasso")
        }
        if( "ridge" %in% unlist(lapply(ls.penalty,"[[","penalty")) ){
            penalty.mg$Vridge <- initVcoef.multigroup(mg, type = "ridge")
        }
        if( "group lasso" %in% unlist(lapply(ls.penalty,"[[","penalty")) ){            
            penalty.mg$Vgroup <- initVcoef.multigroup(mg, type = "group lasso")
        }

        mg$index.numeric <- intersect(manifest(mg$lvm[[1]], lp = FALSE), manifest(x[[1]], lp = FALSE))
        mg$penalty <- penalty.mg

        class(mg) <- c("pmultigroup", class(mg))
    }    
    return(c(list(mg=mg),dots))   
}
# }}}

# {{{ lavaPenalty.post.hook

#' @export
lavaPenalty.post.hook <- function(x){

    
    if(any(c("EPSODE","proximal gradient") %in% x$opt$algorithm )){

        #### update x        
        x$regularizationPath <- x$opt$regPath
        x$opt$regPath <- NULL
        if("multigroupfit" %in% class(x)){
            x$penalty <- x$model$penalty
            x$penaltyNuclear <- x$model$penaltyNuclear            

            x$model$penalty <- NULL
            x$model$penaltyNuclear <- NULL
            
        }else {
            x$penalty <- x$control$penalty
            x$penaltyNuclear <- x$control$penaltyNuclear

            x$control$penalty <- NULL
            x$control$penaltyNuclear <- NULL
        }
        class(x) <- append("plvmfit", class(x))
        
        #### update lambda
        trace <- x$control$trace
        equivariance <- x$control$proxGrad$equivariance
        name.variance <- x$control$proxGrad$name.variance
        fit <- x$control$regPath$fit

        if(is.path(x)){
            ## estimate the best model according to the fit parameter
            if(!is.null(fit)){
                x <- calcLambda(x, fit = fit, trace = trace)                                
            }
      
        }else{
            if(equivariance){
                x$penalty$lambda1.abs <- x$penalty$lambda1
                x$penalty$lambda2.abs <- x$penalty$lambda2
                x$penalty$lambda1 <- x$penalty$lambda1/sum(coef(x)[name.variance])
                x$penalty$lambda2 <- x$penalty$lambda2/sum(coef(x)[name.variance])
            }else{
                x$penalty$lambda1.abs <- x$penalty$lambda1*sum(coef(x)[name.variance])
                x$penalty$lambda2.abs <- x$penalty$lambda2*sum(coef(x)[name.variance])
            }
        }
            
    }
  
    return(x)
  
}

# }}}



#----------------------------------------------------------------------
### Penalty_hook.R ends here
