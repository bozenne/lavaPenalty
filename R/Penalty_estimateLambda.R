#' @title Estimate the regularization parameter
#' 
#' @description Find the optimal regularization parameter according to a fit criterion
#'  
#' @param path the penalized latent variable
#' @param model the non penalized latent variable
#' @param seq_lambda the sequence of penalisation paramaters to be used to define the sub model
#' @param data.fit the data used to fit the model
#' @param data.test the data used to test the model
#' @param warmUp should the new model be initialized with the solution of the previous model (i.e. estimated for a different penalization parameter)
#' @param keep.fit should the penalized LVM be exported
#' @param refit.pLVM should the penalized LVM be fitted on only the non 0 parameters
#' @param fit criterion to decide of the optimal model to retain among the penalized models.
#' @param order the order of the sequence of penalties: according "lambda1" or "lambda1.abs"
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
#'  
#' @export
calcLambda <- function(path, model, 
                       seq_lambda1, data.fit, data.test, 
                       warmUp = FALSE, CI.coef = FALSE,
                       fit = "BIC", order = "lambda1", trace = TRUE, ...){
  
  #### preparation
  if(!missing(path)){
    
    if("lvm" %in% class(model) == FALSE){
      stop("calcLambda: argument \'model\' must be an object of class \"lvm\" \n")
    }
    
    regPath <- getPath(path, getLambda = c("lambda1","lambda1.abs"), getCoef = "coef0", order = order)
    seq_lambda1 <- unlist(lapply(regPath, function(x){attr(x,"lambda1")}))
    seq_lambda1.abs <- unlist(lapply(regPath, function(x){attr(x,"lambda1.abs")})) 
    seq_row <- unlist(lapply(regPath, function(x){attr(x,"row")}))
    seq_coef <- lapply(regPath, function(x){as.character(x)})
    
    penCoef <- path$penCoef
    
  }else{
    
    if("plvm" %in% class(model) == FALSE){
      stop("calcLambda: argument \'model\' must be an object of class \"plvm\" when argument \'path\' is missing \n")
    }
    if(missing(seq_lambda1)){
      stop("calcLambda: argument \'seq_lambda1\' must not be missing when argument \'path\' is missing \n")
    }
    
    path <- list(path = NULL,
                 increasing = all(diff(seq_lambda1)>0),
                 penCoef = model$penalty$name.coef,
                 performance = NULL,
                 optimum = NULL
    )
    class(path) <- "regPath"
    seq_lambda1.abs <- NA
    seq_row <- NA
  }
  
  if(fit %in% c("AIC","BIC","P_error") == FALSE){
    stop("fit must be in AIC BIC P_error \n",
         "proposed value: ",fit,"\n")
  }
  
  if("control" %in% names(list(...))){
    control <- list(...)$control
  }else{
    control <- list()
  }
  if("trace" %in% names(control) == FALSE){
    control$trace <- FALSE
  }
  
  #### initialisation
  n.lambda1 <- length(seq_lambda1)
  n.endogeneous <- length(endogenous(model))
  
  cv.lvm <- rep(NA,n.lambda1)
  seq.criterion <- rep(NA,n.lambda1)
  best.res <- Inf
  best.lambda1 <- NA
  best.subset <- NULL
  best.lambda1.abs <- NULL
  
  fitSave <- NULL
  if(CI.coef){
    pathCI <- list(estimation = NULL, inf = NULL, sup = NULL)
  }
  
  if(trace){pb <- txtProgressBar(min = 0, max = n.lambda1, style = 3)}
  for(iterLambda in 1:n.lambda1){
    
    #### define the variables to include in the model
    if(!missing(path)){
      coef0_lambda <- seq_coef[[iterLambda]]
    }else{
      control1 <- control
      if(warmUp && !is.null(fitSave)){control1$start <- coef(fitSave)}
      fitSave <- estimate(model, data = data.fit, lambda1 = seq_lambda1[iterLambda], control = control1)
      
      path$path <- rbind(path$path, c(lambda1.abs = seq_lambda1[iterLambda],
                                      lambda1 = NA, 
                                      lambda2.abs = NA, 
                                      lambda2 = NA, 
                                      indexChange = NA, 
                                      coef(fitSave)))
      
      coef0_lambda <- coef0(fitSave, tol = 1e-6, penalized = TRUE, value = FALSE)
    }
    
    #### form the reduced model
    model2 <- model
    if(length(coef0_lambda)>0){
      ToKill <- coef0_lambda
      for(iterKill in ToKill){
        model2 <- rmLink.lvm(model2, iterKill)  
      }
    }
    
    #### fit the reduced mode
    fitTempo2 <- lava:::estimate.lvm(model2, data = data.fit, 
                                     control = list(constrain = TRUE, trace = FALSE))
    cv.lvm[iterLambda] <- fitTempo2$opt$convergence==0
    
    if(CI.coef){
      control2 <- control
      control2$start <- coef(fitSave)
      pmodel2 <- penalize(model2, intersect(coef(model2), path$penCoef))
      fitTempo3 <- estimate(pmodel2, data = data.fit, lambda1 = seq_lambda1[iterLambda], control = control2)
      
      pathCI$estimation <- rbind(pathCI$estimation, fitTempo3$coef[,1])
      pathCI$inf <- rbind(pathCI$inf, fitTempo3$coef[,1] - 1.96*fitTempo3$coef[,2])
      pathCI$sup <- rbind(pathCI$sup, fitTempo3$coef[,1] + 1.96*fitTempo3$coef[,2])
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
      best.lambda1 <- seq_lambda1[iterLambda]
      best.subset <- names(coef(fitTempo2))
      best.lvm <- fitTempo2
      
      if(!missing(path)){
        attr(best.lambda1,"row") <- seq_row[iterLambda]
        best.lambda1.abs <- attr(regPath[[as.numeric(attr(best.lambda1,"row"))]],"lambda1.abs")
      }
    }
    
    if(trace){setTxtProgressBar(pb, iterLambda)}
  }
  
  if(trace){close(pb)}
  
  #### export
  path$performance <- data.frame(lambda1 = seq_lambda1,
                                 lambda1.abs = seq_lambda1.abs,
                                 value = seq.criterion,
                                 optimum = (seq.criterion == best.res),
                                 cv = cv.lvm,
                                 row = seq_row)
  
  path$optimum <- list(criterion = fit,
                       value = best.res,
                       coef = coef(best.lvm),
                       lvm = best.lvm,
                       lambda1 = best.lambda1[[1]],
                       lambda1.abs = best.lambda1.abs,
                       row = attr(best.lambda1, "row"))
  
  return(path)
  
}
