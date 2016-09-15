#' @title step of a proximal gradient algorithm
#' @param start initial values for the parameters
#' @param proxOperator proximal operator corresponding to the penalization applied to the log likelihood
#' @param hessian second derivative of the likelihood given by lava. Only used to estimate the step parameter of the algorithm when step = NULL
#' @param gradient first derivative of the likelihood given by lava. 
#' @param objective likelihood given by lava. Used to adjust the step parameter when using backtracking
#' @param lambda1 L1 penalization parameter
#' @param lambda2 L2 penalization parameter
#' @param group.lambda1 group to which each parameter belongs. 0 mean individual lasso otherwise parameters are groupes according to their group.lambda1 value
#' @param step maximun step for the proximal gradient algorithm. 
#' If NULL the step is estimated using the inverse of the maximal eigenvalue of the hessian (in absolute value) and re-estimated at each step
#' Otherwise backtracking is used.
#' @param BT.n number of backtracking steps
#' @param BT.eta multiplicative factor for the step 
#' @param iter.max maximum number of iterations
#' @param abs.tol convergence is the difference in likelihood between two consecutive steps is below this threshold
#' @param rel.tol convergence is the relative difference  in likelihood between two consecutive steps is below this threshold
#' @param fast type of iteration
#' 0 correspond to the ISTA step as described in Bech 2009
#' 1 correspond to the FISTA step as described in Bech 2009
#' 2 correspond to the Monotone APG as described in Li 2015
#' 3 correspond to the Nesterov step as described in Simon 2013
#' @param trace should the convergence diagnostics be displayed at each step
#' 
#' @references 
#' Bech and Teboulle - 2009 A Fast Iterative Shrinkage-Thresholding Algorithm
#' Li 2015 - Accelerated Proximal Gradient Methods for Nonconvex Programming
#' Simon 2013 - A sparse group Lasso
proxGrad <- function(start, proxOperator, method, hessian, gradient, objective,
                     iter.max, trace,
                     abs.tol = lava.options()$proxGrad$abs.tol,
                     rel.tol = lava.options()$proxGrad$rel.tol, 
                     step = lava.options()$proxGrad$step, 
                     BT.n = lava.options()$proxGrad$BT.n, 
                     BT.eta = lava.options()$proxGrad$BT.eta, 
                     force.descent = lava.options()$proxGrad$force.descent,
                     export.iter = lava.options()$proxGrad$export.iter){

  stepMax <- step 
  stepMin <- step*BT.eta^BT.n
  fct_errorLv <- function(e){warning("unable to compute the value of the likelihood - return Inf \n");return(Inf)}
  
  ## initialisation
  x_k <- start 
  
  obj.x_k <- tryCatch(objective(x_k), error = fct_errorLv)
  if(is.na(obj.x_k)){obj.x_k <- Inf}
  grad.x_k <- try(gradient(x_k))
  
  t_k <- t_kp1 <- if(method %in% c("FISTA_Beck")){1}else{NA}
  y_k <- if(method %in% c("FISTA_Beck","FISTA_Vand","mFISTA_Vand")){x_k}else{NA} 
  obj.y_k <- if(method %in% c("FISTA_Beck","FISTA_Vand","mFISTA_Vand")){obj.x_k}else{NA} 
  grad.y_k <- if(method %in% c("FISTA_Beck","FISTA_Vand","mFISTA_Vand")){grad.x_k}else{NA} 

  if("function" %in% class(hessian)){
    maxEigen <- 1/rARPACK::eigs_sym(hessian(x_k),k=1, which = "LM", opts = list(retvec = FALSE))$values
    step <-  abs(maxEigen)
  }
  
  test.cv <- FALSE
  iter <- 0
  iterAll <- 0
  
  if(trace>0){cat("stepBT"," ","iter_back", " ", "max(abs(x_kp1 - x_k))"," ","obj.x_kp1 - obj.x_k","\n")}
  if(export.iter){details.cv <- NULL}
  
    ## loop
  while(test.cv == FALSE && iter <= iter.max){
    iter <- iter + 1 
    
    iter_back <- 0
    diff_back <- 1
    obj.x_kp1 <- +Inf
    
    while( (iter_back < BT.n) && (is.infinite(obj.x_kp1) || diff_back > 0) ){ # Backtracking
      stepBT <- step*BT.eta^iter_back
      iterAll <- iterAll + 1
      
      if(method == "ISTA"){
        res <- ISTA(x_k = x_k, obj.x_k = obj.x_k, grad.x_k = grad.x_k, 
                    proxOperator = proxOperator, step = stepBT)
      }else if(method %in% c("FISTA_Beck","FISTA_Vand","mFISTA_Vand")){
        res <- ISTA(x_k = y_k, obj.x_k = obj.y_k, grad.x_k = grad.y_k, 
                    proxOperator = proxOperator, step = stepBT)
      }
      
      obj.x_kp1 <- tryCatch(objective(res$x_kp1), error = fct_errorLv) 
      if(is.na(obj.x_kp1)){obj.x_kp1 <- Inf}
      
      if(force.descent == TRUE){
        diff_back <- obj.x_kp1 - obj.x_k
      }else{
        diff_back <- obj.x_kp1 - res$Q
      }
      
      iter_back <- iter_back + 1
    
      # cat("obj.x_kp1:",obj.x_kp1," | obj.x_k:",obj.x_k, " | res$Q:",res$Q,"\n")
    }
    
    if(method == "mFISTA_Vand"){
      res$u <- res$x_kp1
      if(obj.x_kp1>obj.x_k){
        res$x_kp1 <- x_k
        obj.x_kp1 <- obj.x_k
        res$cv <- FALSE
      }
    }
    
    # if(obj.x_kp1 > obj.x_k){browser()}
    if(force.descent && obj.x_kp1 > obj.x_k){break}
    
    absDiff <- abs(obj.x_kp1 - obj.x_k) < abs.tol
    relDiff <- abs(obj.x_kp1 - obj.x_k)/abs(obj.x_kp1) < rel.tol
    test.cv <- (absDiff + relDiff > 0)
    if("cv" %in% names(res)){test.cv <- res$cv}
   
    
    #### update
    if(method %in% c("FISTA_Beck","FISTA_Vand","mFISTA_Vand")){
        
      if(method == "FISTA_Beck"){
        t_k <- t_kp1
        t_kp1 <- (1 + sqrt(1 + 4 * t_k^2)) / 2
        y_k <- res$x_kp1 + (t_k-1)/t_kp1 * (res$x_kp1 - x_k) 
      }else if(method == "FISTA_Vand"){
        y_k <- res$x_kp1 + (iter-2)/(iter+1) * (res$x_kp1 - x_k) 
      }else if(method == "mFISTA_Vand"){
        theta_kp1 <- 2/(iter+1)
        v_kp1 <- x_k + 1/theta_kp1 * (res$u - x_k)
        y_k <- (1 - theta_kp1) * res$x_kp1 + theta_kp1 * v_kp1
      }
      
      obj.y_k <- tryCatch(objective(y_k), error = fct_errorLv)
      if(is.na(obj.y_k)){obj.y_k <- Inf}
      grad.y_k <- try(gradient(y_k))
      step <- min(stepMax, stepBT)#min(stepMax, stepBT/sqrt(BT.eta))#
      
    }else{
      
      step <- min(stepMax, stepBT)#min(stepMax, stepBT/sqrt(BT.eta))#
      
    }
    if(trace>0){cat("|",stepBT," ",iter_back, " ", max(abs(res$x_kp1 - x_k))," ",obj.x_kp1 - obj.x_k,"\n")}
    if(export.iter){
      details.cv <- rbind(details.cv,
                          c(iteration = iter, stepBT = stepBT, iter_back = iter_back, adiff_param = max(abs(res$x_kp1 - x_k)), obj = obj.x_kp1, diff_obj = obj.x_kp1 - obj.x_k))
    }
    
    x_k <- res$x_kp1
    obj.x_k <- obj.x_kp1
    grad.x_k <- try(gradient(res$x_kp1))
  }
  if(trace>0){cat("\n")}
  ## export
  message <- if(test.cv){"Sucessful convergence \n"
  }else{
    paste("max absolute/relative difference: ",max(abs(obj.x_kp1 - obj.x_k)),"/",max(abs(obj.x_kp1 - obj.x_k)/abs(obj.x_kp1))," for parameter ",which.max(absDiff),"/",which.max(relDiff),"\n")
  }
  
  return(list(par = x_k,
              step = stepBT,
              convergence = as.numeric(test.cv==FALSE),
              iterations = iter,
              iterationsAll = iterAll,
              evaluations = c("function" = 0, "gradient" = iter),
              message = message,
              details.cv = if(export.iter){details.cv}else{NULL}
  ))
}


ISTA <- function(x_k, obj.x_k, grad.x_k,
                 proxOperator, step){
  
  ## Step
  x_kp1 <- proxOperator(x = x_k - step * grad.x_k, step = step)
  
  ## Upper bound for backtracking
  Q <- Qbound(diff.xy = x_kp1 - x_k, obj.y = obj.x_k, grad.y = grad.x_k, L = 1/step)
  
  return(list(x_kp1 = x_kp1,
              Q = Q))
}

#' @title Estimate an upper bound of obj.x
Qbound <- function(diff.xy, obj.y, grad.y, L){
  
  return(obj.y + crossprod(diff.xy, grad.y) + L/2 * crossprod(diff.xy))
  
}