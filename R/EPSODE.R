# {{{ EPSODE
#' @title Perform the generic Path Algorithm for a LVM
#'
#' @description Perform the generic Path Algorithm for a LVM
#' 
#' @param start the starting value
#' @param objective likelihood given by lava. Used to adjust the step parameter when using backtracking
#' @param gradient first derivative of the likelihood given by lava. 
#' @param hessian second derivative of the likelihood given by lava. Only used to estimate the step parameter of the algorithm when step = NULL
#' @param V matrix that left multiply beta to define the penalization (identity corresponds to a standard lasso penalty)
#' @param lambda2 ridge penalization parameter
#' @param index.penalty2 parameters to which ridge penalization is applied
#' 
#' @param equivariance should the lambda parameter be multiplied with the first variance parameter?
#' @param constrain.variance should the variance parameters be log transformed?
#' @param index.variance the position of the variance parameters in start
#' @param control settings for the EPSODE algorithm. See lava.options.
#'
#' @details Does not work for an unknown variance matrix since the log-likelihood of the regression and variance parameter is not jointly convex
#' (according to ??)
#' 
#' @references 
#' Zhou 2014 - A generic Path Algorithm for Regularized Statistical Estimation


EPSODE <- function(start,
                   objective, gradient, hessian,
                   V,
                   lambda2, index.penalty2,
                   equivariance, constrain.variance, index.variance,
                   control){

    n.coef <- length(start)
    name.coef <- names(start)
    n.penalty <- NCOL(V)
    tol.0 <- control$tol.0
    nstep_max <- control$nstep_max
    
    # {{{ initialize lambda
    if(control$increasing){
        min.0 <- min( abs( ( start %*% V ) / (- gradient(start) %*% V ) ) )
        seq_lambda1 <- 0 
        stepLambda1 <- min.0
        if(is.null(control$stopLambda)){control$stopLambda <- 1e5}
    }else{
        max.grad <- max( abs(gradient(start) %*% V ) )
        seq_lambda1 <-  max.grad * 1.1 # initialisation with the fully penalized solution
        stepLambda1 <- -seq_lambda1
        if(is.null(control$stopLambda)){control$stopLambda <- 0}
    }
    resolution_lambda1 <- control$resolution_lambda1
    # }}}
  
    # {{{ test
    if(length(resolution_lambda1) < 2){
        stop("EPSODE : argument \'resolution_lambda1\' must have length 2 \n",
             "proposed length: ",length(resolution_lambda1),"\n")}
    # }}}
  
    # {{{ constrain variance
    if(equivariance){
        if(constrain.variance){ # on the log scale
            start[index.variance] <- start[index.variance] - start[index.variance[1]]
        }else{ # on the original scale
            start[index.variance] <- start[index.variance]/start[index.variance[1]]
        }
        indexAllCoef <- setdiff(1:n.coef, index.variance[1])
    }else{
        indexAllCoef <- 1:n.coef
    }
    # }}}
  
    # {{{ initialization
    envir <- environment()
    iter <- 1
    test.ncv <- TRUE
  
    ## res
    M.beta <- rbind(start)
    set.penalty <- 1:n.penalty
    setNE <- Matrix::which(M.beta %*% V < -tol.0)
    setPE <- Matrix::which(M.beta %*% V > tol.0)
    setZE <- setdiff(set.penalty, c(setNE, setPE))
    seq_index <- NA
  
    if(control$trace>=1){
        cat("Penalisation path using the EPSODE algorithm \n", sep = "")
        if(equivariance){
            cat(" * fixed coef : \"",paste(names(start)[index.variance[1]], collapse = "\" \""),"\" \n", sep = "")
        }
        if(control$trace==1){pb <- utils::txtProgressBar(min = 0, max = n.penalty, style = 3)}
    }
    # }}}

    # {{{ main loop
    while(iter < nstep_max && test.ncv>0){

        ### current parameters
        iterLambda1 <- tail(seq_lambda1,1)
        if(iter > 1 && iterLambda1 == seq_lambda1[iter-1]){
            iterLambda1 <- iterLambda1 + stepLambda1 * resolution_lambda1[2]
        }
        iterBeta <- M.beta[NROW(M.beta),]
    
        ### estimate uz and Uz
        uz <- rep(0, n.coef)
        if(length(setNE)>0){uz <- uz - Matrix::rowSums(V[,setNE,drop = FALSE])}
        if(length(setPE)>0){uz <- uz + Matrix::rowSums(V[,setPE,drop = FALSE])}
        Uz <- as.matrix(Matrix::t(V[,setZE,drop = FALSE]))#[setZE,indexAllCoef,drop = FALSE]

        if(length(setZE)>0){ # in case of ill conditionned problem
            B <- pracma::nullspace(Uz[,indexAllCoef,drop = FALSE])
            iUz <- MASS::ginv(Uz[,indexAllCoef,drop = FALSE]) #  solve(Uz %*% t(Uz)) %*% Uz #  or (Uz_pen t(Uz_pen))^-1 Uz_pen
        }else{
            B <- NULL
            iUz <- NULL
        }
    
        ### Solve ODE 
        lambda.ode <- seq_len(500)
        cv.ode <- NULL
        bridge.ode <- rbind(c(iter = 0, lambda = iterLambda1,  step = 0, iterBeta))

        res.error <- try(deSolve::ode(y = iterBeta,
                                      times = lambda.ode,
                                      func = EPSODE_odeBeta, method = control$ode.method,
                                      parm = list(hessian = hessian, 
                                                  setNE = setNE, setZE = setZE, setPE = setPE, 
                                                  lambda2 = lambda2, index.penalty2 = index.penalty2, indexAllCoef = indexAllCoef, 
                                                  V = V, uz = uz, Uz = Uz, iUz = iUz, B = B,
                                                  resolution = resolution_lambda1, lars = control$lars, reversible = control$reversible, increasing = control$increasing, stopLambda = control$stopLambda,
                                                  envir = envir)
                                      ), silent = TRUE)
        
        ## remove initialization row and duplicated rows
        index.duplicated <- which(!duplicated(bridge.ode[-1,c(-1,-3), drop = FALSE]))+1
        bridge.ode <- bridge.ode[index.duplicated,,drop=FALSE]
        
        ## update
        if(control$exportAllPath && length(bridge.ode[,"iter"])>2){ ## export all the points of the regularization path
            seq_index <- c(seq_index, rep(NA, nrow(bridge.ode)-2) )       
            seq_lambda1 <- c(seq_lambda1, 
                             unname(bridge.ode[c(-1,-nrow(bridge.ode)),"lambda",drop = FALSE]))
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
        }else if(cv.ode$cv["cv.lambda0"]){
            seq_index <- c(seq_index, NA)
        }

        newLambda1 <- unname(bridge.ode[nrow(bridge.ode),"lambda"])
        newBeta <- unname(bridge.ode[nrow(bridge.ode),name.coef])
        if(length(setZE)>0){newBeta[setZE] <- 0}
    
        M.beta <- rbind(M.beta, newBeta)
        seq_lambda1 <- c(seq_lambda1, newLambda1)
        iter <- iter + 1
    
        ## cv
        if(is.na(newLambda1)){break}
        if(control$increasing){
            test.ncv <- (length(setNE) > 0 || length(setPE) > 0 )
            if(!is.null(control$stopLambda) && newLambda1>=control$stopLambda){test.ncv <- -1}
            if(!is.null(control$stopParam) && length(setZE)>=control$stopParam){test.ncv <- -1}
        }else {
            test.ncv <- newLambda1>0
            #test.ncv <- length(setZE)>0
            #if(newLambda1==0){test.ncv <- -1 ;  seq_index <- c(seq_index, NA) ;}
            if(newLambda1<=control$stopLambda && control$stopLambda>0){test.ncv <- -1}
            if(!is.null(control$stopParam) && (length(setNE)+length(setPE))>=control$stopParam){test.ncv <- -1}
        }
        if(control$trace>1){
            cat("iteration ",iter-1,": lambda=",newLambda1,"\n") ; print(bridge.ode[nrow(bridge.ode),name.coef])
        }else if(control$trace==1){
            utils::setTxtProgressBar(pb, value = if(control$increasing){length(setZE)}else{length(setNE)+length(setPE)})
        }
    
    }
    # }}}

    if(control$trace==1){close(pb)}
    
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

    dt <- as.data.table(cbind(index = 1:NROW(M.beta),
                              lambda1.abs = if(equivariance == 0){NA}else{seq_lambda1}, 
                              lambda1 = if(equivariance == 0){seq_lambda1}else{NA}, 
                              lambda2.abs = if(!is.null(lambda2)){mean(lambda2)}else{NA}, 
                              lambda2 = NA,
                              indexChange = unname(seq_index),
                              M.beta))
    
    return(list(par = NULL,
                message = message,
                convergence =  as.numeric(test.ncv),
                iterations = iter,
                algorithm = "EPSODE",
                path = dt))
}

# }}}

# {{{ EPSODE_odeBeta
EPSODE_odeBeta <- function(t, y, ls.args){

    bridge <- get("bridge.ode", envir = ls.args$envir)
    lambda1 <- sum(bridge[tail(which(bridge[,"iter"]==(t-1)),1),c("lambda","step")])
    equivariance1 <- as.vector(y %*% ls.args$V)

    # {{{ test lambda exceed stopLambda
    if((ls.args$increasing==FALSE && lambda1 <= ls.args$stopLambda) || (ls.args$increasing && lambda1 >= ls.args$stopLambda)){
   
        assign("cv.ode", 
               value = list(param = y, lambda = lambda1, index = NA, cv = c(cv.sign = FALSE, cv.constrain = FALSE, cv.lambda0 = TRUE, s = NA)), 
               envir = ls.args$envir)
        assign("bridge.ode",
               value =  rbind(bridge,c(t, lambda = ls.args$stopLambda, 0, y)),
               envir = ls.args$envir)
        stop("EPSODE_odeBeta: lambda1 = stopLambda \n",
             "end of the path \n")
    }
    # }}}
    
    # {{{ check constrains: sign change
    if(ls.args$reversible || ls.args$increasing){
        index <- NULL
        if(length(ls.args$setNE)>0){index <- c(index, ls.args$setNE[which(equivariance1[ls.args$setNE] > 0)])} ## any negative constrain that becomes positive: stop algorithm
        if(length(ls.args$setPE)>0){index <- c(index, ls.args$setPE[which(equivariance1[ls.args$setPE] < 0)])} ## any positive constrain that becomes negative: stop algorithm
    
        if(length(index) == 1){
            assign("cv.ode",
                   value = list(param = y, lambda = lambda1, index = index, cv = c(cv.sign = TRUE, cv.constrain = FALSE, cv.lambda0 = FALSE, s = NA)),
                   envir = ls.args$envir)
            assign("bridge.ode",
                   value =  rbind(bridge,c(t, lambda1, 0, y)),
                   envir = ls.args$envir)
            stop("cv \n")
        }else if(length(index) > 1){
            assign("cv.ode",
                   value = list(param = y, lambda = lambda1, index = NA, cv = c(cv.sign = length(index), cv.constrain = FALSE, cv.lambda0 = FALSE, s = NA)),
                   envir = ls.args$envir)
            assign("bridge.ode",
                   value =  rbind(bridge,c(t, NA, NA, y)),
                   envir = ls.args$envir)
            stop("EPSODE_odeBeta: multiple events \n",
                 "increase resolution \n")
        }
    }
    # }}}
    
    # {{{ Hessian and gradient
    Hfull <- ls.args$hessian(y)
    H <- Hfull[ls.args$indexAllCoef,ls.args$indexAllCoef, drop = FALSE]
    attr(H, "grad")  <- attr(Hfull, "grad")[ls.args$indexAllCoef, drop = FALSE]
    if(!is.null(ls.args$lambda2)){
        attr(H, "grad")[ls.args$index.penalty2] <- attr(H, "grad")[ls.args$index.penalty2] + ls.args$lambda2 * y[ls.args$index.penalty2, drop = FALSE]
        H[ls.args$index.penalty2,ls.args$index.penalty2] <- H[ls.args$index.penalty2,ls.args$index.penalty2] + diag(ls.args$lambda2)
    }  
    G <- attr(H, "grad")
    uz <- ls.args$uz[ls.args$indexAllCoef, drop = FALSE]
    Uz <- ls.args$Uz[,ls.args$indexAllCoef, drop = FALSE]
    # }}}

    # {{{ estimate Q and P
    if(length(ls.args$setZE) == 0){
        R <- NULL
        Q <- NULL
        P <- solve(H)
    }else{
        H_m1 <- try(solve(H), silent = TRUE)
        if(is.matrix(H_m1)){ # 
            # all coef
            R <- - solve(Uz %*% H_m1 %*% t(Uz))
            Q <- H_m1 %*% t(Uz) %*% (-R)            
            P <- H_m1 - Q %*% Uz %*% H_m1 
         
        }else{ # singular H matrix
            if(all(H<1e-12)){
                assign("cv.ode", 
                       value = list(param = y, lambda = lambda1, index = NA, cv = c(cv.sign = FALSE, cv.constrain = FALSE, cv.lambda0 = FALSE, s = NA)), 
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
        # }}}

        # {{{ check constrains: subgradient
        if(ls.args$reversible || ls.args$increasing==FALSE){
            s <- - t(Q) %*% ( (1 / lambda1) * G + uz) # - t(Q) %*% ( (1 / nextKnot) * G + uz)
            if(t > 1 && any(abs(s) - 1 > -abs(ls.args$resolution[2])) ){
                index <- which.max(abs(s))
                assign("cv.ode",
                       value =  list(param = y, lambda = lambda1, index = ls.args$setZE[index], cv = c(cv.sign = FALSE, cv.constrain = TRUE, s = s[index])),
                       envir = ls.args$envir)
                assign("bridge.ode",
                       value =  rbind(bridge,c(t, lambda1, 0, y)),
                       envir = ls.args$envir)
                stop("cv \n")
            }
        }
    }
    # }}}

    # if(any(uz)>0) browser()
    # {{{ Compute differential of the ODE
    Puz <- rep(0, length(y))
    Puz[ls.args$indexAllCoef] <- P %*% uz

    #uz["Y~X2"] <- 0
    # P2 <- P
    # rownames(P2) <- names(y)
    # colnames(P2) <- names(y)
    # P2
    
    ## remove noise due to numeric approximations
    ## Puz <- round(Puz,12)
 
    # }}}

    # {{{ Compute the new lambda
    ## control the evolution of each parameter
    if(ls.args$lars || (t==1 && length(c(ls.args$setNE,ls.args$setPE))==0)){ # no limit
        # when all parameter are shrinked to 0, the second member of the ODE is constant with lambda so linear interpolation is ok
        absDiff.deltaLambda <- Inf
        
    }else{ # the change should be ls.args$resolution[1] between two steps

        coef.n0 <- setdiff(ls.args$indexAllCoef,ls.args$setZE)
        absDiff.deltaLambda <- min(ls.args$resolution[1]/abs(Puz[coef.n0]))
    }
   
    ## find next node
    if(ls.args$increasing){ # when a coefficient become 0

        possible.deltaLambda <- as.numeric(y/Puz) %*% ls.args$V[,c(ls.args$setNE,ls.args$setPE),drop=FALSE]
        li.deltaLambda <- min(possible.deltaLambda[possible.deltaLambda>=0])
        deltaLambda <- max(min(li.deltaLambda,absDiff.deltaLambda),ls.args$resolution[2])
        
    }else{ # when the constrains won't be respected
        if(!is.null(Q)){
            lPlus <- (- t(Q) %*% G)/(1 + t(Q) %*% uz)
            lMinus <- (- t(Q) %*% G)/(-1 + t(Q) %*% uz)
            lAll <- c(lPlus[lPlus<lambda1],lMinus[lMinus<lambda1])
            li.deltaLambda <- lAll[which.min(lambda1-lAll)]-lambda1          
           
        }else{
            li.deltaLambda <- -lambda1
        }
         deltaLambda <- min(max(li.deltaLambda,-absDiff.deltaLambda),-ls.args$resolution[2])
    }
    # }}}

    ## export
    assign("bridge.ode",
           value =  rbind(bridge,c(t, lambda1, deltaLambda, y)),
           envir = ls.args$envir)
     
    return(list(-Puz*deltaLambda))
    #    return(list(-Puz*normTempo))
}
# }}}

