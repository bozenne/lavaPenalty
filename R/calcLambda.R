# {{{ doc
#' @title Estimate the regularization parameter
#' @name calcLambda
#' @description Find the optimal regularization parameter according to a fit criterion
#'  
#' @param x a plvmfit object
#' @param seq_lambda1 the sequence of penalisation paramaters to be used to define the sub model
#' @param data.fit the data used to fit the model
#' @param data.test the data used to test the model
#' @param warmUp should the new model be initialized with the solution of the previous model (i.e. estimated for a different penalization parameter)
#' @param fit criterion to decide of the optimal model to retain among the penalized models.
#' @param trace shoud the execution of the function be traced
#' @param ... additional arguments - e.g. control argument for estimate.lvm
#' 
#' @examples 
#' 
#' set.seed(10)
#' n <- 300
#' formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:4), collapse = "+")))
#' mSim <- lvm(formula.lvm)
#' df.data <- sim(mSim,n)
#' 
#' lvm.model <- lvm(formula.lvm)
#' plvm.model <- penalize(lvm.model)
#'
#' res <- estimate(plvm.model, data = df.data, increasing = FALSE, regularizationPath = TRUE)
#' 
#' perf <- calcLambda(res$regularizationPath, res$x, data.fit = df.data)
#' perf

#' @rdname calcLambda
#' @export
`calcLambda` <- function(x, ...) UseMethod("calcLambda")
# }}}

# {{{ calcLambda.plvm
#' @rdname calcLambda
#' @export
calcLambda.plvmfit <- function(x, 
                               seq_lambda1,
                               data.fit = x$data$model.frame,
                               data.test = x$data$model.frame, 
                               warmUp = lava.options()$calcLambda$warmUp, 
                               fit = lava.options()$calcLambda$fit,
                               trace = TRUE,
                               ...){

    # {{{ test
    if(is.path(x)==FALSE){
        stop("missing regularization path in object \n")
    }
    if(fit %in% c("AIC","BIC","P_error") == FALSE){
        stop("fit must be in AIC BIC P_error \n",
             "proposed value: ",fit,"\n")
    }
    
    # }}}

    # {{{ extract path
    constrain <- penalty(x, type = "Vlasso")$Vlasso
        
    regPath <- getPath(x, path.constrain = TRUE, only.breakpoints=TRUE)
    seq_lambda1.abs <- regPath$lambda1.abs
    seq_lambda1 <- regPath$lambda1
    n.lambda1 <- length(seq_lambda1.abs)
    # }}}
    
    
    # {{{ initialization
    if("control" %in% names(list(...))){
        control <- list(...)$control
    }else{
        control <- list()
    }
    if("trace" %in% names(control) == FALSE){
        control$trace <- FALSE
    }
      
    n.knot <- NROW(regPath)
    n.constrain <- NCOL(penalty(x, type = "Vlasso")$Vlasso)
    n.coef <- length(coef(x))

    cv.lvm <- rep(NA,n.knot)
    seq.criterion <- rep(NA,n.knot)
    best.res <- Inf
    best.lambda1 <- NA
    best.lambda1.abs <- NULL
  
    pathCI <- data.table(expand.grid(coef = names(coef(x)), knot = 1:n.knot))
    pathCI[, estimate := as.numeric(NA)]
    pathCI[, lower := as.numeric(NA)]
    pathCI[, upper := as.numeric(NA)]
    pathCI[, p.value := as.numeric(NA)]
    pathCI[, lambda1 := regPath$lambda1[knot]]
    pathCI[, lambda2 := regPath$lambda2[knot]]
    pathCI[, lambda1.abs := regPath$lambda1.abs[knot]]
    pathCI[, lambda2.abs := regPath$lambda2.abs[knot]]
    setkeyv(pathCI, c("coef","knot"))

    store.resCoef <- setNames(rep(NA, n.coef), names(coef(x)))
    
    ## attribute a name to each coefficient
    model0 <- x$model
    coefWithPenalty <- names(Matrix::which(Matrix::rowSums(abs(constrain))>0))
    index.coefWithPenalty <- setNames(match(coefWithPenalty, rownames(constrain)),coefWithPenalty)

    for(iCoef in coefWithPenalty){ # iCoef <- coefWithPenalty[1]
        regression(model0, as.formula(iCoef)) <- as.formula(paste0("~beta", index.coefWithPenalty[iCoef]))
    }
    # }}}
    
    if(trace){
        cat("Estimation of the optimum lambda according to the ",fit," criterion \n")
        pb <- txtProgressBar(min = 0, max = n.lambda1, style = 3)
    }

    
    for(iKnot in 1:n.knot){ # 
        
        ## define the constrains
        ils.constrain <- regPath$constrain0[[iKnot]]        
        
        # {{{ form the reduced model
        Imodel0 <- model0
        class(Imodel0) <- setdiff(class(model0),"plvm")
        if(!is.null(ils.constrain)){      
            for(iCon in ils.constrain){ # iCon <- 1
                iConstrain <- constrain[,iCon]
                iName <- names(which(iConstrain!=0))
                iCoef <- iConstrain[iConstrain!=0]
                
                if(length(iName)==1){
                    regression(Imodel0, as.formula(iName)) <- 0
                }else{
                    iF <- as.formula(paste0(index.coefWithPenalty[iName[1]],"~",
                                            index.coefWithPenalty[iName[-1]]))
                    constrain(Imodel0, iF) <- function(x){}
                }
            }
        }
        # }}}
        
        # {{{ fit the reduced model and extract gof        
        fitTempo2 <- estimate(Imodel0, data = data.fit, 
                              control = control)
        ## estimates
        tableTempo2 <- summary(fitTempo2)$coef
        
        store.resCoef[] <- NA
        store.resCoef[rownames(tableTempo2)] <- tableTempo2[,"Estimate"]
        pathCI[knot == iKnot, estimate := store.resCoef]
        
        store.resCoef[] <- NA
        store.resCoef[rownames(tableTempo2)] <- tableTempo2[,"Estimate"] + qnorm(0.025)*tableTempo2[,"Std. Error"]
        pathCI[knot == iKnot, lower := store.resCoef]
        
        store.resCoef[] <- NA
        store.resCoef[rownames(tableTempo2)] <- tableTempo2[,"Estimate"] + qnorm(0.975)*tableTempo2[,"Std. Error"]
        pathCI[knot == iKnot, upper := store.resCoef]
        
        ## gof criteria
        cv.lvm[iKnot] <- fitTempo2$opt$convergence==0
        if(fit %in% c("AIC","BIC")){
            seq.criterion[iKnot] <- gof(fitTempo2)[[fit]]
        }else{
            predY <- as.data.frame(predict(fitTempo2,
                                           data = data.test[,exogenous(fitTempo2), drop = FALSE]))
      
            seq.criterion[iKnot] <- 1/(n.endogeneous*NROW(data.test)) * sum((data.test[,endogenous(fitTempo2)] - predY)^2)
        }

        ## look for the best LVM
        if(cv.lvm[iKnot] && best.res>seq.criterion[iKnot]){
            best.res <- seq.criterion[iKnot]
            best.lvm <- fitTempo2
            best.lambda1 <- seq_lambda1[iKnot]
            best.lambda1.abs <-seq_lambda1.abs[iKnot]
        }
        # }}}
        
        if(trace){setTxtProgressBar(pb, iKnot)}
    }

    if(trace){close(pb)}
    
    
    #### export
    best.lvm$control <- x$control
    best.lvm$call <- x$call
    
    best.lvm$regularizationPath <- x$regularizationPath
    best.lvm$penalty <- x$penalty
    best.lvm$penalty$lambda1 <- best.lambda1
    best.lvm$penalty$lambda1.abs <- best.lambda1.abs
    best.lvm$nuclearPenalty <- x$nuclearPenalty
    
    best.lvm$regularizationPath$path[match(index,regPath$index), (fit) := seq.criterion]
    best.lvm$regularizationPath$path[match(index,regPath$index), cv := cv.lvm]
    best.lvm$regularizationPath$path[match(index,regPath$index), optimum := seq.criterion==best.res]
    best.lvm$regularizationPath$criterion <- fit

    class(best.lvm) <- class(x)
    return(best.lvm)
    
}
# }}}

