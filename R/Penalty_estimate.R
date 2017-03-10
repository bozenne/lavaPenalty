
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
estimate.plvm <- function(x, data, 
                          lambda1 = NULL, lambda2 = NULL, lambdaG = NULL, lambdaN = NULL, adaptive = FALSE, 
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

    # {{{ prepare data (scaling, convert character to dummy variables)
    if(control$trace>0){cat("Scale and center dataset \n")}
    resData <- prepareData.plvm(x, data = data)
    x <- resData$lvm
    data <- resData$data
    # }}}  
    # {{{ prepare penalty
    ### update penalty (for new variables)
    x <- updatePenalty.plvm(x, data = data)

    ### update the penalty with potential user input
    if(!is.null(lambda1)){
        penalty(x, type = "lambda1") <- lambda1
    }
    if(!is.null(lambda2)){    
        penalty(x, type = "lambda2") <- lambda2
    }
    if(!is.null(lambdaG)){    
        penalty(x, type = "lambdaG") <- lambdaG
    }
    if(!is.null(lambdaN)){    
        penalty(x, type = "lambdaN", nuclear = TRUE) <- lambdaN
    }
    if(!is.null(adaptive)){    
        penalty(x, type = "adaptive") <- adaptive
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

# {{{ updatePenalty.penaltyL12
#' @title Prepare the data for the estimate function
#'  
#' @description Scale the data, update the penalty term according to the presence of factors. 
#' 
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param method.center function used to center the data
#' @param method.scale function used to scale the data
#' 
updatePenalty.plvm <- function(x, data){

    #### elastic net
    link.penalty <- penalty(x, type = "link")
    if(any(link.penalty %in% coef(x) == FALSE)){
        missing.links <- link.penalty[link.penalty %in% coef(x) == FALSE]

        for(iLink in missing.links){ # iLink <- missing.links[1]
            resInit  <- lava.reduce::initVar_link(iLink)
            if(class(data[[resInit$var2]]) == "factor"){
                cancelPenalty(x) <- iLink
            }

        }
    }
    group.penalty <- penalty(x, type = "group")
    
    
    #### convert categorical variables to dummy variables
    resC2D <- lava:::categorical2dummy(x, data)
  
    index.numeric <- intersect(manifest(x, lp = FALSE), manifest(resC2D$x, lp = FALSE))
    indexOld.factor <- setdiff(manifest(x, lp = FALSE),  manifest(resC2D$x, lp = FALSE))
    indexNew.factor <- setdiff(manifest(resC2D$x, lp = FALSE), manifest(x, lp = FALSE))
  
    test.factor <- length(indexNew.factor)>0
  
    #### update the penalty according to the dummy variables
    name.Newlinks <- coef(x)
    browser()
    # find the old coefficients corresponding to the factor variables
    OldLinksPenalty.factors <- setdiff(link.penalty, name.Newlinks)
    
    if(test.factor && length(OldLinksPenalty.factors)>0){
        
        #     find the new coefficient to penalize
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
              scale = value.scale,
              center = value.center,
              lvm = x))
}
# }}}
