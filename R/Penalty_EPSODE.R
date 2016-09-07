#' @title Perform the generic Path Algorithm for a LVM
#' 
#' 
#' @param beta_lambda0 values of the parameters when the penalty parameter is 0.
#' @param beta_lambdaMax values of the parameters when the penalty parameter is infinite
#' @param objective likelihood given by lava. Used to adjust the step parameter when using backtracking
#' @param gradient first derivative of the likelihood given by lava. 
#' @param hessian second derivative of the likelihood given by lava. Only used to estimate the step parameter of the algorithm when step = NULL
#' @param V matrix that left multiply beta to define the penalization (identity corresponds to a standard lasso penalty)
#' @param lambda2 L2 penalization parameter
#' 
#' @param indexPenalty position of the penalised coefficients in beta
#' @param indexNuisance index of the nuisance parameter to be treated as a constant
#' 
#' @param resolution_lambda1 the first value is the maximum relative difference in parameter between two steps. 
#' If this lead to a too small step, the second value is used as the minimum change in penalization parameter between two steps.
#' @param increasing direction of the path
#' @param stopLambda if not null, stop the path when the penalty parameter reach this value
#' @param stopParam if not null, stop the path when the number of 0 (if increasing = TRUE) or non 0 (if increasing = FALSE) parameters has reached this value.
#' @param nstep_max the maximum number of iterations
#' @param ode.method the type of method to use to solve the ode (see the documentation of deSolve:::ode)
#' @param control additional options to be passed to the proximal algorithm
#' @param reversible should the algorithm allow a 0 parameter to become non-0 (when increasing = TRUE) or non-0 parameter to become 0 (when increasing = FALSE)
#' @param tol.0 tolerance for classify a parameter from beta_lambdaMax in the set of 0 parameters
#' @param exportAllPath export all the regularization path (and not only the breakpoints)
#' @param trace should the function be traced
#' 
#' @references 
#' Zhou 2014 - A generic Path Algorithm for Regularized Statistical Estimation


EPSODE <- function(beta_lambda0, beta_lambdaMax, objective, gradient, hessian, V, lambda2, 
                   indexPenalty, indexNuisance, 
                   resolution_lambda1, increasing, stopLambda, stopParam,
                   nstep_max = min(length(beta)*50,1e4), 
                   ode.method = "euler", control, reversible, tol.0 = 1e-8, exportAllPath, trace){
  
  #### preparation
  if(increasing){
    beta <- beta_lambda0
  }else{
    beta <- beta_lambdaMax
  }
  
  ## lambda
  res <- initLambda_EPSODE(increasing = increasing,
                           gradient = gradient, beta = beta_lambdaMax, indexPenalty = indexPenalty, indexNuisance = indexNuisance)
  seq_lambda1 <- res$seq_lambda
  stepLambda1 <- res$stepLambda
  if(increasing == FALSE){resolution_lambda1 <- -resolution_lambda1}
  if(length(resolution_lambda1) < 2){
    stop("EPSODE : argument \'resolution_lambda1\' must have length 2 \n",
         "proposed length: ",length(resolution_lambda1),"\n")}
  
  
  n.coef <- length(beta)
  lambda2_save <- lambda2
  lambda2 <- rep(0, n.coef)
  lambda2[indexPenalty] <- lambda2_save 
  envir <- environment()
  
  ## constrain 
  if(length(indexNuisance) > 0){
    res <- initSigmaConstrain(beta, constrain = control$constrain, indexNuisance = indexNuisance)
    beta <- res$start
    indexAllCoef <- res$indexAllCoef
  }else{
    constrain <- NULL
    indexAllCoef <- 1:n.coef
  }
  
  #### initialization
  iter <- 1
  test.ncv <- TRUE
  
  ## res
  M.beta <- rbind(beta)
  setNE <- intersect(which(V %*% beta < -tol.0),  indexPenalty)
  setZE <- intersect(which(abs(V %*% beta) < tol.0), indexPenalty)
  setPE <- intersect(which(V %*% beta > tol.0), indexPenalty)
  seq_index <- NA
  
  if(trace>=0){
    cat("Penalisation path using the EPSODE algorithm \n", sep = "")
    if(length(indexNuisance) > 0){
      cat(" * fixed coef : \"",paste(names(beta)[indexNuisance], collapse = "\" \""),"\" \n", sep = "")
    }
    if(trace==0){pb <- utils::txtProgressBar(min = 0, max = length(indexPenalty), initial = 0, style = 3)}
  }
  
  #### main loop
  while(iter < nstep_max && test.ncv>0){
    
    ## current parameters
    iterLambda1 <- tail(seq_lambda1,1)
    if(iter > 1 && iterLambda1 == seq_lambda1[iter-1]){
      iterLambda1 <- iterLambda1 + stepLambda1*resolution_lambda1[2]
    }
    iterBeta <- M.beta[nrow(M.beta),]
    
    #### estimate uz and Uz
    uz <- rep(0, n.coef)
    if(length(setNE)>0){uz <- uz - colSums(V[setNE,,drop = FALSE])}
    if(length(setPE)>0){uz <- uz  + colSums(V[setPE,,drop = FALSE])}
    Uz <- V[setZE,indexAllCoef,drop = FALSE]
    
    if(length(Uz)>0){ # in case of ill conditionned problem
      B <- pracma::nullspace(Uz)
      iUz <- solve(Uz %*% t(Uz)) %*% Uz # MASS::ginv(Uz_pen) # or (Uz_pen t(Uz_pen))^-1 Uz_pen
    }else{
      B <- NULL
      iUz <- NULL
    }
    
    ## Solve ODE 
    lambda.ode <- seq_len(1000)
    cv.ode <- NULL
    bridge.ode <- rbind(c(iter = 0, step = 0, lambda = iterLambda1, iterBeta))
    
    # H1 <- hessian(iterBeta)
    # H2 <- hessianGaussianO(iterBeta)
    # Hdiff <- H1 - H2
    # attr(Hdiff,"grad") <- attr(H1,"grad") - attr(H2,"grad")
    # print(Hdiff)
    res.error <- try(deSolve::ode(y = iterBeta,
                                  times = lambda.ode,
                                  func = EPSODE_odeBeta, method = ode.method,
                                  parm = list(hessian = hessian, setNE = setNE, setZE = setZE, setPE = setPE,
                                              lambda2 = lambda2, indexPenalty = indexPenalty, indexAllCoef = indexAllCoef, indexNuisance = indexNuisance,
                                              uz = uz, Uz = Uz, iUz = iUz, B = B,
                                              resolution = resolution_lambda1, reversible = reversible, envir = envir)
    ), silent = TRUE)
    
     ## second chance in case of multiple events
    if(!is.null(cv.ode) && cv.ode$cv["cv.sign"]>1){
      bridge.odeS <- bridge.ode
      iterBeta2 <- bridge.ode[max(1,nrow(bridge.odeS)-10),-(1:3)]
      lambda.ode2 <- unname(bridge.ode[max(1,nrow(bridge.odeS)-10),3])
      
      lambda.ode <- seq_len(10000)
      cv.ode <- NULL
      bridge.ode <- rbind(c(iter = 0, step = 0, lambda = lambda.ode2, iterBeta2))
      
      
      res.error <- try(deSolve::ode(y = iterBeta2,
                                    times = lambda.ode,
                                    func = EPSODE_odeBeta, method = ode.method,
                                    parm = list(hessian = hessian, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE,
                                                lambda2 = lambda2, indexPenalty = indexPenalty, indexAllCoef = indexAllCoef,
                                                uz = uz, Uz = Uz, iUz = iUz, B = B,
                                                resolution = resolution_lambda1/10, reversible = reversible, envir = envir)
      ), silent = TRUE)
      
    }
    #  cat("\n iteration ",iter,"\n")

    ## update 
    if(exportAllPath && length(bridge.ode[,"iter"])>2){ ## export all the points of the regularization path
      seq_index <- c(seq_index, rep(NA, nrow(bridge.ode)-2) )  
     
      seq_lambda1 <- c(seq_lambda1, 
                       bridge.ode[c(-1,-nrow(bridge.ode)),"lambda",drop = FALSE])
      M.beta <- rbind(M.beta, 
                      bridge.ode[c(-1,-nrow(bridge.ode)),-(1:3),drop = FALSE])
    }
    
    if(is.null(cv.ode) || cv.ode$cv["cv.sign"]>1){
      if(all(class(res.error) != "try-error")){ ## not enought interations: continue
        seq_index <- c(seq_index, NA)  
      }else{
        stop(res.error[1])
      }
    }else if(cv.ode$cv["cv.sign"]==1){
      setNE <- setdiff(setNE,  cv.ode$index)
      setPE <- setdiff(setPE, cv.ode$index)
      setZE <- union(setZE, cv.ode$index)
      seq_index <- c(seq_index, cv.ode$index)
    }else if(cv.ode$cv["cv.constrain"]){
      setZE <- setdiff(setZE, cv.ode$index)
      if(cv.ode$cv["s"]>0){
        setPE <- union(setPE, cv.ode$index) 
      }else{
        setNE <- union(setNE, cv.ode$index)  
      }
      seq_index <- c(seq_index, cv.ode$index)
    }
    
    newLambda1 <- unname(bridge.ode[nrow(bridge.ode),3])
    newBeta <- unname(bridge.ode[nrow(bridge.ode),-(1:3)])
    if(length(setZE)>0){newBeta[setZE] <- 0}
    
    M.beta <- rbind(M.beta, newBeta)
    seq_lambda1 <- c(seq_lambda1, newLambda1)
    iter <- iter + 1
    
    ## cv
    if(is.na(newLambda1)){break}
    if(stepLambda1 > 0){
      test.ncv <- (length(setNE) > 0 || length(setPE) > 0 )
      if(!is.null(stopLambda) && newLambda1>=stopLambda){test.ncv <- -1}
      if(!is.null(stopParam) && length(setZE)>=stopParam){test.ncv <- -1}
    }else {
      test.ncv <- length(setZE)>0
      if(newLambda1==0){test.ncv <- -1 ;  seq_index <- c(seq_index, NA) ;}
      if(!is.null(stopLambda) && newLambda1<=stopLambda){test.ncv <- -1}
      if(!is.null(stopParam) && (length(setNE)+length(setPE))>=stopParam){test.ncv <- -1}
    }
    if(trace>=1){
      cat("iteration ",iter-1,": lambda=",newLambda1,"\n") ; print(bridge.ode[nrow(bridge.ode),-(1:3)])
    }else if(trace==0){
      utils::setTxtProgressBar(pb, value = if(increasing){length(setZE)}else{length(setNE)+length(setPE)})
    }
    
  }
  if(trace==0){close(pb)}
  
  #### post treatment
  if(increasing == FALSE && !is.null(beta_lambda0) && (0 %in% seq_lambda1 == FALSE) && test.ncv==0){
    M.beta <- rbind(M.beta,
                    unname(beta_lambda0))
    seq_lambda1 <- c(seq_lambda1, 0)
    seq_index <- c(seq_index, NA)
  }
  
  #### export
  if(test.ncv == FALSE){
    message <- "Sucessful convergence \n"
  }else if(iter >= nstep_max){
    message <- "Maximum number of steps reached \n"
  }else if(is.na(newLambda1) || newLambda1 < 0){
    message <- "Invalid penalization parameter \n"
  }
  
  seq_lambda1 <- unname(seq_lambda1)
  rownames(M.beta) <- NULL
  df <- as.data.frame(cbind(lambda1.abs = if(length(indexNuisance) == 0){NA}else{seq_lambda1}, 
                            lambda1 = if(length(indexNuisance) == 0){seq_lambda1}else{NA}, 
                            lambda2.abs = lambda2_save, 
                            lambda2 = NA,
                            indexChange = unname(seq_index),
                            M.beta))
   return(list(message = message,
               path = df))
}




EPSODE_odeBeta <- function(t, y, ls.args){
  
  bridge <- get("bridge.ode", envir = ls.args$envir)
  lambda <- bridge[tail(which(bridge[,1]==(t-1)),1),3]
  
  #### test lambda <= 0
  if(lambda < 0 || ls.args$resolution[2] < 0 && lambda == 0){
    
    assign("cv.ode", 
           value = list(param = y, lambda = lambda, index = NA, cv = c(cv.sign = FALSE, cv.constrain = FALSE, s = NA)), 
           envir = ls.args$envir)
    assign("bridge.ode",
           value =  rbind(bridge,c(t, NA, 0, y)),
           envir = ls.args$envir)
    stop("EPSODE_odeBeta: lambda = 0 \n",
         "end of the path \n")
  }
  
  #### test knot
  if(ls.args$reversible || ls.args$resolution[2]>0){ #### to be removed, here to avoid instabilities
    index <- NULL
    if(length(ls.args$setNE)>0){index <- c(index, ls.args$setNE[which(y[ls.args$setNE] > 0)])} ## any negative parameter that becomes positive: stop algorithm
    if(length(ls.args$setPE)>0){index <- c(index, ls.args$setPE[which(y[ls.args$setPE] < 0)])} ## any positive parameter that becomes negative: stop algorithm
    
    if(length(index) == 1){
      assign("cv.ode",
             value = list(param = y, lambda = lambda, index = index, cv = c(cv.sign = TRUE, cv.constrain = FALSE, s = NA)),
             envir = ls.args$envir)
      assign("bridge.ode",
             value =  rbind(bridge,c(t, NA, lambda, y)),
             envir = ls.args$envir)
      stop("cv \n")
    }else if(length(index) > 1){
      assign("cv.ode",
             value = list(param = y, lambda = lambda, index = NA, cv = c(cv.sign = length(index), cv.constrain = FALSE, s = NA)),
             envir = ls.args$envir)
      assign("bridge.ode",
             value =  rbind(bridge,c(t, NA, NA, y)),
             envir = ls.args$envir)
      stop("EPSODE_odeBeta: multiple events \n",
           "increase resolution \n")
    }
  }
  
  #### Hessian and gradient
  Hfull <- ls.args$hessian(y)
  H <- Hfull[ls.args$indexAllCoef,ls.args$indexAllCoef, drop = FALSE]
  attr(H, "grad")  <- attr(Hfull, "grad")[ls.args$indexAllCoef, drop = FALSE]
  if(any(ls.args$lambda2>0)){
    attr(H, "grad") <- attr(H, "grad") + ls.args$lambda2[ls.args$indexAllCoef, drop = FALSE] * y[ls.args$indexAllCoef, drop = FALSE]
    H[] <- H[] + diag(ls.args$lambda2[ls.args$indexAllCoef, drop = FALSE])
  }
  
  G <- attr(H, "grad")[ls.args$indexAllCoef, drop = FALSE]
  uz <- ls.args$uz[ls.args$indexAllCoef, drop = FALSE]
  Uz <- ls.args$Uz[,ls.args$indexAllCoef, drop = FALSE]
  
  #### estimate Q and P
  if(length(ls.args$setZE) == 0){
    R <- NULL
    Q <- NULL
    P <- solve(H)
  }else{
    H_m1 <- try(solve(H), silent = TRUE)
    if(is.matrix(H_m1)){ # 
      # all coef
      H_m1 <- solve(H)
      R <- solve(Uz %*% H_m1 %*% t(Uz))
      Q <- H_m1 %*% t(Uz) %*% R
      P <- H_m1 - Q %*% Uz %*% H_m1 
         
    }else{ # singular H matrix
      if(all(H<1e-12)){
        assign("cv.ode", 
               value = list(param = y, lambda = lambda, index = NA, cv = c(cv.sign = FALSE, cv.constrain = FALSE, s = NA)), 
               envir = ls.args$envir)
        assign("bridge.ode",
               value =  rbind(bridge,c(t, NA, NA, y)),
               envir = ls.args$envir)
        stop("EPSODE_odeBeta: all values in the hessian are below 1e-12 \n",
             "try using a larger step \n")
      }
      BHB <-  t(ls.args$B) %*% H %*% ls.args$B
      P <- ls.args$B %*% solve(BHB) %*% t(ls.args$B)
      Q <- t(ls.args$iUz[,ls.args$indexAllCoef, drop = FALSE])
      
    }
    
    ## check constrains
    if(ls.args$reversible || ls.args$resolution[2]<0){ #### to be removed, here to avoid instabilities
      s <- - t(Q) %*% ( (1 / lambda) * G + uz) # - t(Q) %*% ( (1 / nextKnot) * G + uz)
      if(t > 1 && any( abs(s) > 1)){
        index <- which.max(abs(s))
        assign("cv.ode",
               value =  list(param = y, lambda = lambda, index = ls.args$setZE[index], cv = c(cv.sign = FALSE, cv.constrain = TRUE, s = s[index])),
               envir = ls.args$envir)
        assign("bridge.ode",
               value =  rbind(bridge,c(t, NA, lambda, y)),
               envir = ls.args$envir)
        stop("cv \n")
      }
    }
  }
  
  ## export
  Puz <- rep(0, length(y))
  if(length(ls.args$setZE)==0){
    coef.n0 <- ls.args$indexAllCoef
    Puz[ls.args$indexAllCoef] <- P %*% uz
  }else{
    coef.n0 <- setdiff(ls.args$indexAllCoef,ls.args$setZE)
    Puz[setdiff(ls.args$indexAllCoef,ls.args$setZE)] <- P[-ls.args$setZE,-ls.args$setZE, drop = FALSE] %*% ls.args$uz[setdiff(ls.args$indexAllCoef,ls.args$setZE), drop = FALSE]  
  }
  
  if(ls.args$resolution[1]>0){
    # lAll <- c(y[ls.args$setNE] / Puz[ls.args$setNE], y[ls.args$setPE] / Puz[ls.args$setPE])
    # if(any(lAll>0)){
    #   ## linear interpolation
    #   normTempo <- max(min(lAll[lAll>0]*ls.args$resolution[1]), ls.args$resolution[2])
    # }else{
    ## relative difference
    rdiff.max <- max(abs(Puz[coef.n0]/y[coef.n0]))
    normTempo <- max(ls.args$resolution[1]/rdiff.max,ls.args$resolution[2])  
    # }
    
  }else{
   ## linear interpolation
    lPlus <- (- t(Q) %*% G)/(1 + t(Q) %*% uz)
    lMinus <- (- t(Q) %*% G)/(-1 + t(Q) %*% uz)
    lAll <- c(lPlus[lPlus<lambda],lMinus[lMinus<lambda])
    nextKnot <- lAll[which.min(lambda-lAll)]
    
    if(t==1 && length(c(ls.args$setNE,ls.args$setPE))==0){ ## when all parameter are shrinked to 0, the second member of the ODE is constant with lambda so linear interpolation is ok
      sign <- c(-1,1)[which.min(c(min(abs(lMinus-nextKnot)), min(abs(lPlus-nextKnot))))]
      if(sign<0){
        index <- ls.args$setZE[which.min(abs(lMinus-nextKnot))]
      }else{
        index <- ls.args$setZE[which.min(abs(lPlus-nextKnot))]
      }
      
      assign("cv.ode",
             value =  list(param = y, lambda = nextKnot,
                           index = index,
                           cv = c(cv.sign = FALSE, cv.constrain = TRUE, s = sign)),
             envir = ls.args$envir)
      assign("bridge.ode",
             value =  rbind(bridge,c(t, NA, nextKnot, y)),
             envir = ls.args$envir)
      stop("cv \n")
      
    }else{
      normTempo <- max(-lambda,min( (lambda-nextKnot)*ls.args$resolution[1], ls.args$resolution[2])) ## max(-lambda) to avoid negative lambda
    }
    
  }
  
  assign("bridge.ode",
         value =  rbind(bridge,c(t, normTempo, lambda+normTempo, y)),
         envir = ls.args$envir)
  return(list(-Puz*normTempo))
}


initLambda_EPSODE <- function(increasing, gradient, beta, indexPenalty, indexNuisance){
  
  if(increasing){
    seq_lambda <- 0 
    
    stepLambda <- max( abs(-gradient(beta) )[indexPenalty] )
    if(length(indexNuisance) > 0){stepLambda <- stepLambda * beta[indexNuisance[1]]}
    
    
  }else{
    
    seq_lambda <-  max( abs(-gradient(beta) )[indexPenalty] ) * 1.1 # initialisation with the fully penalized solution
    if(length(indexNuisance) > 0){seq_lambda <- seq_lambda * beta[indexNuisance[1]]}
    
    stepLambda <- -seq_lambda
    
  }
  
  return(list(seq_lambda = seq_lambda,
              stepLambda = stepLambda))
}


initSigmaConstrain <- function(start, constrain, indexNuisance){
  
  if(constrain){
    start[indexNuisance] <- start[indexNuisance] - start[indexNuisance[1]]
    constrain <- setNames(0, names(start)[indexNuisance[1]])
  }else{
    start[indexNuisance] <- start[indexNuisance]/start[indexNuisance[1]]
    constrain <- setNames(1, names(start)[indexNuisance[1]])
  }
  indexAllCoef <- setdiff(1:length(start), indexNuisance[1])
  
  return(list(start = start,
              constrain = constrain,
              indexAllCoef = indexAllCoef)
  )
}