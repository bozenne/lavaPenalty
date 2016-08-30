calcLambda <- function(object, model, seq_lambda, data.fit, data.test, 
                       warmUp = FALSE, keep.fit = FALSE, refit.pLVM = TRUE,
                       fit = "BIC", order = "lambda1", trace = TRUE, ...){
 
   # if no test set then CV
 if(fit %in% c("AIC","BIC","P_error") == FALSE){
   stop("fit must be in AIC BIC P_error \n")
 }
  if(all(c("plvmfit", "lvm") %in% class(object) == FALSE)){
    stop("model must be either a plvmfit or a lvm object \n")
  }
  
  if("control" %in% names(list(...))){
    control <- list(...)$control
  }else{
    control <- list()
  }
  if("trace" %in% names(control) == FALSE){
    control$trace <- FALSE
  }
  
  if("plvmfit" %in% class(object)){
    regPath <- getPath(object, getLambda = order, getCoef = "coef0", order = order)
    seq_lambda <- unlist(lapply(regPath, function(x){attr(x,order)}))
    seq_row <- unlist(lapply(regPath, function(x){attr(x,"row")}))
    seq_coef <- lapply(regPath, function(x){as.character(x)})
    
    regPath2 <- list(pLVMred = NULL)
  }else{
    if(missing(model)){model <- object}
    ls.model <- NULL
    ls.modelRed <- NULL
    regPath2 <- list(pLVM = NULL, pLVMred = NULL)
  }
  
  ##
  n.lambda <- length(seq_lambda)
  n.endogeneous <- length(endogenous(model))
 
  cv.lvm <- rep(NA,n.lambda)
  seq.criterion <- rep(NA,n.lambda)
  best.res <- Inf
  best.lambda <- NA
  best.subset <- NULL
  
  fitTempo1 <- NULL
  
  if(trace){pb <- txtProgressBar(min = 0, max = n.lambda, style = 3)}
  for(iterLambda in 1:n.lambda){
    
    #### define the variables to include in the model
    if("plvmfit" %in% class(object) == FALSE){
      control1 <- control
      if(warmUp && !is.null(fitTempo1)){control1$start <- coef(fitTempo1)}
      fitTempo1 <- estimate(model, data = data.fit, lambda1 = seq_lambda[iterLambda], control = control1)
      regPath2$pLVM <- rbind(regPath2$pLVM, c(lambda = seq_lambda[iterLambda], 
                                              cv = fitTempo1$opt$convergence==0, coef(fitTempo1)))
      if(keep.fit){
        ls.model <- c(ls.model, list(fitTempo1))
      }
      coef0_lambda <- coef0(fitTempo1, tol = 1e-6, penalized = TRUE, value = FALSE)
    }else{
      coef0_lambda <- seq_coef[[iterLambda]]
    }
    
    #### form the reduced model
    model2 <- model
    if(length(coef0_lambda)>0){
      ToKill <- coef0_lambda
      for(iterKill in ToKill){
        model2 <- rmLink.lvm(model2, iterKill)  
      }
    }
  
    #### fit the reduced model
    fitTempo2 <- lava:::estimate.lvm(model2, data = data.fit, 
                                     control = list(constrain = TRUE, trace = FALSE, start = coef(fitTempo1)))
    cv.lvm[iterLambda] <- fitTempo2$opt$convergence==0
    
    if(refit.pLVM && "plvmfit" %in% class(object) == FALSE){
      control2 <- control
      control2$start <- coef(fitTempo1)
      pmodel2 <- penalize(model2, intersect(coef(model2),model$penalty$name.coef))
      fitTempo3 <- estimate(pmodel2, data = data.fit, lambda1 = seq_lambda[iterLambda], control = control2)
      
      newRow <- regPath2$pLVM[iterLambda,]
      newRow[names(coef(fitTempo3))] <- coef(fitTempo3)
      regPath2$pLVMred <- rbind(regPath2$pLVMred, newRow)
      
      newRow.inf <- regPath2$pLVM[iterLambda,]
      newRow.inf[names(coef(fitTempo3))] <- fitTempo3$coef[,1] - 1.96*fitTempo3$coef[,2]
      regPath2$pLVMred.inf <- rbind(regPath2$pLVMred, newRow.inf)
      
      newRow.sup <- regPath2$pLVM[iterLambda,]
      newRow.sup[names(coef(fitTempo3))] <- fitTempo3$coef[,1] + 1.96*fitTempo3$coef[,2]
      regPath2$pLVMred.sup <- rbind(regPath2$pLVMred, newRow.sup)

      if(keep.fit){
        ls.modelRed <- c(ls.modelRed, list(fitTempo3))
      }
    }
    
    #### gof criteria
    if(fit %in% c("AIC","BIC")){
      seq.criterion[iterLambda] <- gof(fitTempo2)[[fit]]
    }else{
      predY <- as.data.frame(predict(fitTempo2,
                                     data = data.test[,exogenous(fitTempo2), drop = FALSE]))
      
      seq.criterion[iterLambda] <- 1/(n.endogeneous*NROW(data.test)) * sum((data.test[,endogenous(fitTempo2)] - predY)^2)
    }
    
    #### storage
    if(cv.lvm[iterLambda] && best.res>seq.criterion[iterLambda]){
      best.res <- seq.criterion[iterLambda]
      best.lambda <- seq_lambda[iterLambda]
      best.subset <- names(coef(fitTempo2))
      best.lvm <- fitTempo2
      if("plvmfit" %in% class(object)){attr(best.lambda,"row") <- seq_row[iterLambda]}
    }
    
    if(trace){setTxtProgressBar(pb, iterLambda)}
  }
  
  if(trace){close(pb)}
  if("plvmfit" %in% class(object)){
    best.lvm$penalty <- object$penalty 
    best.lvm$penalty$lambda1 <- seq_lambda 
    attr(seq.criterion, "criterion") <- fit
    best.lvm$penalty$performance <- seq.criterion
    best.lvm$penalty$convergence <- cv.lvm
    best.lvm$penalty$lambda1.best <- best.lambda
    best.lvm$regularizationPath <- object$regularizationPath
    class(best.lvm) <- append("plvmfit", class(best.lvm))
    return(best.lvm)
  }else{
    return(list(criterion = seq.criterion,
                lambda = best.lambda,
                subset = best.subset,
                lvm = best.lvm,
                cv = cv.lvm,
                regPath = regPath2,
                plvm.fit = ls.model,
                plvmRed.fit = ls.modelRed))
  }
}
