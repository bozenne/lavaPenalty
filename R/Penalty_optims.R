
# {{{ optim.proxGrad
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
optim.proxGrad <- function(start, objective, gradient, hessian, control, ...){

    ## only keep the relevant elements in control
    control.proxGrad <- control$proxGrad
    penalty <- control$penalty
    penaltyNuclear <- control$penaltyNuclear
    name.coef <- names(start)
    
    # {{{ update according to the coefficients present in start
    initL12 <- initialize.penaltyL12(penalty, name.coef = name.coef,
                                     regularizationPath = FALSE)
    initNuclear <- initializer.penaltyNuclear(penaltyNuclear, name.coef = name.coef)        
    if(!identical(initNuclear$lambdaN,0)){
        hessian <- NULL
    }
    index.variance <- match(control.proxGrad$name.variance, name.coef)
    # }}}
    
    # {{{ Proximal gradient
    ## generate proximal operator
    resInit <- initializeOperator(lambda1 = initL12$lambda1, index.penalty1 = initL12$index.penalty1,
                                  lambda2 = initL12$lambda2, index.penalty2 = initL12$index.penalty2,
                                  lambdaG = initL12$lambdaG, index.penaltyG = initL12$index.penaltyG,
                                  lambdaN = initNuclear$lambdaN, index.penaltyN = initNuclear$index.penaltyN, 
                                  nrow = initNuclear$nrow, ncol = initNuclear$ncol,
                                  constrain.lambda = control.proxGrad$constrain.lambda,
                                  constrain.variance = control.proxGrad$constrain.variance,
                                  index.variance = index.variance)
    # resInit$objectivePenalty(start)
    # resInit$proxOperator(start, step = 0.1)
    
    ## run descent
    res <- proxGrad(start = start, proxOperator = resInit$proxOperator,
                    hessian = hessian, gradient = gradient, objective = objective, 
                    control = control.proxGrad)
    # }}}
    
    # {{{ Adaptive proximal gradient
    if(control.proxGrad$adaptive){
        ## generate proximal operator
        resInit <- initializeOperator(lambda1 = initL12$lambda1/abs(res$par),
                                      lambda2 = initL12$lambda2,
                                      lambdaG = initL12$lambdaG, index.penaltyG = initL12$index.penaltyG,
                                      lambdaN = initNuclear$lambdaN, index.penaltyN = initNuclear$index.penaltyN, 
                                      nrow = initNuclear$nrow, ncol = initNuclear$ncol,
                                      constrain.lambda = control$proxGrad$fixSigma,
                                      constrain.variance = control$proxGrad$constrain,
                                      index.variance = control$proxGrad$name.variance)

        ## run descent    
        res <- proxGrad(start = res$par, proxOperator = resInit$proxOperator,
                        hessian = hessian, gradient = gradient, objective = objective, 
                        control = control.proxGrad)
    }

    # }}}

    ## update objective with penalty
    res$objectivePenalty <- resInit$objectivePenalty

    ## export
    return(res)
}
# }}}

# {{{ optim.proxGradPath
optim.proxGradPath <- function(start, objective, gradient, hessian, control, ...){

    ## extract
    trace <- control$trace
    grid <- as.data.frame(control$regPath$grid)
    warmUp <- control$regPath$warmUp
    control$regPath <- NULL
    control$proxGrad$trace <- FALSE
    
    ## init
    n.grid <- NROW(grid)
    name.grid <- names(grid)

    name.coef <- names(start)
    n.coef <- length(start)
    res.Path <- data.table(matrix(as.numeric(NA), nrow = n.grid, ncol = n.coef+7))
    names(res.Path) <- c("index","lambda1.abs","lambda1","lambda2.abs","lambda2","indexChange",name.coef,"cv")

    if(trace){
        cat("Penalisation path using the proxGrad algorithm \n", sep = "")
        pb <- txtProgressBar(min = 1, max = n.grid)
    }
    
    for(iLambda in 1:n.grid){ # iLambda <- 1
        ## update penalty
        for(iPen in 1:NCOL(grid)){ # iPen <- 1
            penalty(control$penalty, type = name.grid[[iPen]], add = FALSE) <- grid[iLambda,iPen]
            res.Path[iLambda,(name.grid[[iPen]]) := grid[iLambda,iPen]]
        }

        ## update start
        if(warmUp && iLambda>1){
            iStart <- res$par
        }else{
            iStart <- start
        }
        
        ## estimation
        res <- optim.proxGrad(start, objective, gradient, hessian, control, ...)

        ## store results        
        vec.res <- c(iLambda,res$par,as.numeric(res$convergence==0))
        vec.names <- c("index",name.coef,"cv")
        res.Path[iLambda, (vec.names) := as.list(vec.res)]
        if(trace){setTxtProgressBar(pb, iLambda)}
    }

    if(trace){close(pb)}

    ## export
    regPath <- list(path = res.Path,
                    increasing = TRUE,
                    criterion = NULL)
    class(regPath) <- "regPath"

    if(all(res.Path$cv==1)){
        message <- "Sucessful convergence \n"
    }else{
        message <- "Convergence issues \n"
    }
    
    res <- list(par = start,#setNames(rep(NA, length(start)), names(start)),
                convergence = 0,
                iterations = 0,
                algorithm = "proximal gradient",
                message = message, 
                regPath = regPath,
                objective = NA)

}
# }}}

# {{{ optim.regEPSODE

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
optim.EPSODE <- function(start, objective, gradient, hessian, control, ...){

    name.coef <- names(start)
    penalty <- control$penalty
    control.regPath <- control$regPath

    # {{{ update according to the coefficients present in start
    initL12 <- initialize.penaltyL12(penalty, name.coef = name.coef,
                                     regularizationPath = TRUE)
    index.variance <- match(control.regPath$name.variance, name.coef)
    # }}}

    # {{{ EPSODE
    resEPSODE <- EPSODE(start = start,
                        objective = objective, gradient = gradient, hessian = hessian, 
                        V = initL12$Vlasso,
                        lambda2 = initL12$lambda2, index.penalty2 = initL12$index.penalty2,
                        constrain.lambda = control.regPath$constrain.lambda,
                        constrain.variance = control.regPath$constrain.variance,
                        index.variance = index.variance,
                        control = control.regPath)
    # }}}

    # {{{ estimation of the nuisance parameter and update lambda/lambda.abs
    if(control.regPath$constrain.lambda){
        res <- optim.Nuisance.plvm(x = control$regPath$envir$x,
                                   data = control$regPath$envir$data,
                                   path = resEPSODE$path[,name.coef,with=FALSE],
                                   index.nuisance = index.variance,
                                   control = control,
                                   ...)   

        if(length(index.variance)==1){
            sigma2 <- resEPSODE$path[[control.regPath$name.variance]]
            if(control.regPath$constrain.variance){sigma2 <- exp(sigma2)}
            resEPSODE$path[, lambda1.abs := lambda1.abs*sigma2]
            resEPSODE$path[, lambda2.abs := lambda2.abs*sigma2]
        }

        resEPSODE$path[,control.regPath$name.variance := lapply(1:NCOL(res),function(c){res[,c]})]

        if(length(index.variance)==1){            
            if(control.regPath$constrain.variance){res <- exp(res)}            
            resEPSODE$path[, lambda1 := lambda1.abs/res]
            resEPSODE$path[, lambda2 := lambda2.abs/res]
        }
    }else{
        if(length(index.variance)==1){
            resEPSODE$path[, lambda1.abs := lambda1*resEPSODE$path[[control.regPath$name.variance]]]
            resEPSODE$path[, lambda2.abs := lambda2*resEPSODE$path[[control.regPath$name.variance]]]
        }
    }
    # }}}
    
    # {{{ conversion to regPath object
    regPath <- list(path = resEPSODE$path,
                    increasing = control.regPath$increasing,
                    criterion = NULL)
    class(regPath) <- "regPath"
    # }}}

    res <- list(par = start,#setNames(rep(NA, length(start)), names(start)),
                convergence = resEPSODE$convergence,
                iterations = resEPSODE$iterations,
                algorithm = resEPSODE$algorithm,
                message = resEPSODE$message, 
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
                                path, index.nuisance,
                                control, ...){

    ## clean control
    trace <- control$trace
    control$trace <- FALSE
    control$regPath <- NULL
    control$proxGrad <- NULL
    control$penalty <- NULL
    control$penaltyNuclear <- NULL
    control$start <- NULL
    control$method <- NULL
    
    ## prepare index coef
    n.coef <- NCOL(path)
    n.knots <- NROW(path)
    n.nuisance <- length(index.nuisance)
    names.coef <- names(path)
    index.constrain <-  setdiff(1:n.coef,c(1, index.nuisance))
    
    ## constrain LVM model
    class(x) <- setdiff(class(x),"plvm")

    if("lvm.reduced" %in% class(x)){
        stop("optim.Nuisance needs to be updated to deal with lvm.reduce object [TODO] \n")
       # to be updated for removing lp when we need to fix links
        if(names.coef[iter_p] %in% lp(xConstrain, type = "link")){
            x <- cancelLP(x, link = names.coef[iP])
            x <- addLink(x, var1 = names.coef[iP], covariance = FALSE)
        }
    }

    ##  estimate nuisance for each node
    if(trace>=1){cat("Estimation of the nuisance parameter ")}

    M.nuisance <- matrix(NA, nrow = n.knots, ncol = n.nuisance)

    for(iKnot in 1:n.knots){ # iKnot <- 1
        xConstrain <- x
        for(iP in index.constrain){ # iP <- 1
            xConstrain <- setLink(xConstrain, var1 = names.coef[iP], value = path[iKnot][[iP]])    
        }# path[iKnot,]       
        elvm <- estimate(xConstrain, data = data, control = control)
        if( abs(path[iKnot,1]-coef(elvm)[1] ) > 1e-3 ){
            warning("possible inaccuracy when estimating the nuisance parameter \n")
        }
        M.nuisance[iKnot,] <- coef(elvm)[names.coef[index.nuisance]]
    }

    if(trace>=1){cat("- done \n")}

    ## export
    return(M.nuisance)  
}

# }}}



