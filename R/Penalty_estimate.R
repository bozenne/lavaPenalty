
# {{{ estimate

# {{{ doc
#' @title Penalized lvm model
#' @name estimate
#' @aliases estimate estimate.plvm
#' @description Estimate a penalized lvm model
#'
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param lambda1 lasso penalization parameter
#' @param lambda2 ridge penalization parameter
#' @param lambdaG group lasso penalization parameter
#' @param lambdaN nuclear norm penalization parameter
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
#' R/examples/EX_Penalty_estimate.R
# }}}

# {{{ estimate.plvm
#' @rdname estimate
#' @export
estimate.plvm <- function(x, data, lambda1, lambda2, lambdaG, lambdaN, adaptive = FALSE, 
                          regularizationPath = FALSE, fit = lava.options()$calcLambda$fit,
                          control = list(), control.proxGrad = list(), control.EPSODE = list(), 
                          constrain.lambda = FALSE, ...) {

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
                          constrain.variance = control$constrain,
                          constrain.lambda = constrain.lambda,
                          name.variance = coefVar(x, value = TRUE)
                          )
    if (length(controlUser.proxGrad) > 0) { # update with user input
        control$proxGrad[names(controlUser.proxGrad)] <- controlUser.proxGrad
    }

    ## regularization path
    if(regularizationPath){   
    
        control$regPath <- c(lava::lava.options()$EPSODE,
                             start.lambda0 = NULL,
                             start.lambdaMax = NULL,
                             fit = fit)

        # if linear regression then use LARS algorithm
        test.regression <- (length(endogenous(x)) == 1) && (length(latent(x)) == 0)
        if(test.regression && all( c("resolution_lambda1","increasing") %in% names(controlUser.regPath) == FALSE)){
            control$regPath$resolution_lambda1 <- c(1,1e-10)
            control$proxGrad$constrain.lambda <- TRUE
            control$proxGrad$increasing <- FALSE
        }
        
        if (length(controlUser.regPath) > 0) { # update with user input
            control$regPath[names(controlUser.regPath)] <- controlUser.regPath
        }
    }else{
        control$regPath <- NULL
    }
    # }}}


    # {{{ prepare penalty
    ### update penalty (for factor variables)
    x <- updatePenaltyFactor.plvm(x, data = data,
                                  trace = control$trace)
    table.penalty <- penalty(x, nuclear = FALSE)

    ### update the penalty with potential user input
    if("lasso" %in% table.penalty$penalty){
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
    resData <- prepareData.plvm(x, data = data)
    x <- resData$lvm
    data <- resData$data
    
    
    if(any(table.penalty$link %in% coef(x) == FALSE)){
        stop("unknown penalized links \n",
             "links: \"",paste(table.penalty$link[table.penalty$link %in% coef(x) == FALSE], collapse = "\" \""),"\"\n")
    }    
    # }}}
    
   
    
    # {{{ optimisation
    res <- lava_estimate.lvm(x = x, data = data,
                             method = if(regularizationPath == 0){"optim.regLL"}else{"optim.regPath"},
                             control = control, quick = FALSE, index = TRUE, ...)
    # }}}
    
    ## export
    return(res)
}
# }}}

# }}}

# {{{ updatePenalty.plvm
#' @title Update the lasso penalty to group lasso penalty in presence of categorical variables
#'  
#' @description Update the lasso penalty to group lasso penalty in presence of categorical variables
#' 
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param trace should the user be told that some penalties have been updated ?
#'
#' @examples
#'
#' m <- lvm(Y~X1+X2+X3+X4)
#' pm <- penalize(m)
#' updatePenaltyFactor.plvm(pm, data = sim(m, 1), trace = TRUE)
#' 
#' mCAT <- m
#' categorical(mCAT, K = 3, labels = letters[1:3]) <- "X1"
#' categorical(mCAT, K = 2, labels = letters[1:2]) <- "X2"
#'
#' pm <- penalize(m)
#' updatePenaltyFactor.plvm(pm, data = sim(mCAT, 1), trace = TRUE)
#'
#' pmCAT <- penalize(mCAT)
#' updatePenaltyFactor.plvm(pmCAT, data = sim(mCAT, 1), trace = TRUE)
#' 
updatePenaltyFactor.plvm <- function(x, data, trace){

    table.penalty <- penalty(x, nuclear = FALSE) 
    var.penalty <- unique(table.penalty$exogenous)        
    test.factor <- sapply(var.penalty, function(v){is.factor(data[[v]])})
    
    ## search for links that should be moved from lasso to group lasso
    if(any(test.factor)){
       
        var.factor <- unique(var.penalty[test.factor])
        table.penalty.factor <- table.penalty[table.penalty$exogenous %in% var.factor & penalty == "lasso"]

        ## remove lasso link
        cancelPenalty(x, lasso = TRUE, ridge = FALSE, group = FALSE) <- table.penalty.factor[["link"]]

        ## update the lvm with the categorical variables
        for(f in var.factor){ # f <- "X1"
            if(trace>0){cat("convert the lasso penalty on the categorical variable ",f," to a group penalty \n")}
            categorical(x, labels = levels(as.factor(data[[f]]))) <- as.formula(paste0("~",f))
        }
            
        ## update the penalty according to the categorical variables
        resInit <- initGroup.lvm(x, links = table.penalty.factor$link)
        # same as `penalize<-.lvm`
        newV <- initVcoef.lvm(x,
                              link = resInit[type == "categorical"][["link"]],
                              group = resInit[type == "categorical"][["group"]])
        penalty(x, type = "Vgroup", add = TRUE) <- newV
    }

    #### export
    return(x)
}
# }}}

# {{{ updatePenaltyVariable.penaltyL12
#' @title Remove penalties corresponding to reference links
#'  
#' @description Remove penalties corresponding to reference links
#' 
#' @param x a penalized lvm model
#' @param name.coef the names of the coefficients used in the optimisation routine
updatePenaltyVariable.penaltyL12 <- function(x, name.coef){
    n.coef <- length(name.coef)
    table.penalty <- penalty(x)

    ## check 
    if(any(table.penalty$link %in% name.coef == FALSE)){
        rm.penalty <- setdiff(table.penalty$link, name.coef)
        warning("initPenalty: some penalty will not be applied because the corresponding parameter is used as a reference \n",
                "non-applied penalty: ",paste(rm.penalty, collapse = " "),"\n")
        cancelPenalty(x, rm.lasso = TRUE, rm.ridge = TRUE, rm.group = TRUE, extraParameter = NULL) <- rm.penalty
     }

    return(x)
}

# }}}



# {{{ initializer.penaltyNuclear
#' @rdname initializer
initializer.penaltyNuclear <- function(x, name.coef){
    n.coef <- length(name.coef)
    name.allPenalty <- penalty(x, type = "link")
    lambdaN <- penalty(x, type = "lambdaN")
    if(length(name.allPenalty)==0){
        return(list(lambdaN = 0, ncol = NULL, nrow = NULL, index.penaltyN = NULL))
    }

    n.nuclear <- length(name.allPenalty)
    vec.lambdaN <- setNames(rep(0, n.coef), name.coef)
    vec.ncol <- penalty(x, type = "ncol")
    vec.nrow <- penalty(x, type = "nrow")
    index.penaltyN <- vector(mode = "list", length= n.nuclear)
    
    for(iNuclear in 1:n.nuclear){

        if(any(name.allPenalty[[iNuclear]] %in% name.coef == FALSE)){
            stop("initPenalty: some penalty will not be applied because the corresponding parameter is used as a reference \n",
                 "non-applied penalty: ",paste(setdiff(name.allPenalty[[iNuclear]], name.coef), collapse = " "),"\n")
        }
        vec.lambdaN[name.allPenalty[[iNuclear]]] <- lambdaN[iNuclear]
        index.penaltyN[[iNuclear]] <- setNames(match(name.allPenalty[[iNuclear]],
                                                     name.coef), name.allPenalty[[iNuclear]])

    
    }
        
    return(list(lambdaN = vec.lambdaN, ncol = vec.ncol, nrow = vec.nrow,
                index.penaltyN = index.penaltyN))
}
# }}}


# {{{ prepareData.plvm
#' @title Scale the dataset corresponding to a LVM object.
#'  
#' @description Scale the data and update the type of the variables in the LVM object according to the variables in the dataset.
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
  
    #### export
    return(list(data = data,
                conversion.factor = conversion.factor,
                scale = value.scale,
                center = value.center,
                lvm = x))
}
# }}}
