#### objective / gradient / hessian
objectiveMC <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  n <- length(Y)
  Xint <- cbind(1,X)
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  mc <-  t(epsilon) %*% epsilon
  
  return(as.numeric(mc))
}

gradientMC <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  n <- length(Y)
  Xint <- cbind(1,X)
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  gradient_sigma2 <- 0
  gradient_beta <- - t(Xint) %*% epsilon  
  
  return(as.numeric(c(gradient_beta,gradient_sigma2)))
}

hessianMC <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  G <- gradientMC(coef = coef, Y = Y, X = X)
  
  n <- length(Y)
  Xint <- cbind(1,X)
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  hessian_beta <- t(Xint) %*% Xint
  
  hessian_sigma2 <- 0
  hessian_sigma2FDbeta <- rep(0, length(beta))
  
  H <- cbind( rbind(hessian_beta,t(hessian_sigma2FDbeta)),
              c(hessian_sigma2FDbeta,hessian_sigma2)
  )
  attr(H,"grad") <- G
  
  return(-H)
}


coef2.penalized <- function(x, iter_lambda, name.response){
  
  if(is.list(x)){
    
    if(!missing(iter_lambda)){
      x <- x[[iter_lambda]]
    }else{
      res <- lapply(x, function(model){
        c(model@lambda1,
          model@lambda2,
          model@unpenalized, 
          model@penalized, 
          model@nuisance$sigma2)
      })
      Mres <- matrix(unlist(res), nrow = length(res), byrow = TRUE)
      colnames(Mres) <- c("lambda1","lambda2",names(x[[1]]@unpenalized),names(x[[1]]@penalized),"sigma2")
      return(Mres)
    }
    
  }
  
    coef <- c(x@unpenalized, 
              x@penalized, 
              x@nuisance$sigma2)
    n.coef <- length(coef)
    if(missing(name.response)){name.response <- all.vars(x@formula$penalized)[1]}
    names(coef)[1] <- name.response
    names(coef)[2:(n.coef-1)] <- paste(name.response, names(coef)[2:(n.coef-1)], sep = lava.options()$symbols[1])
    names(coef)[length(names(coef))] <- paste(name.response, name.response, sep =  lava.options()$symbols[2])
    return(coef)
}

validPath.lvm <- function(x, data, validMean = TRUE, validCov = TRUE, control = list(), ...){
  
  Mall <- getPath(x, lambda = c("lambda1","lambda2"), only.breakpoints = TRUE)
  Mcoef <- Mall[,-(1:2),drop = FALSE]
  seqLambda1 <- Mall[,"lambda1"]
  seqLambda2 <- Mall[,"lambda2"]
  seqLambda2[is.na(seqLambda2)] <- 0
  n.Lambda <- length(seqLambda1)
  
  if("trace" %in% names(control) == FALSE){
    control$trace <- FALSE
  }
  
  McoefGS <- NULL
  for(iterL in 1:n.Lambda){
    McoefGS <- rbind(McoefGS, coef(estimate(x$x, data = data, control = control,
                                        lambda1 = seqLambda1[iterL], lambda2 = seqLambda2[iterL],
                                        ...)))
  }
  
  return(list(proxAlgo = McoefGS,
              EPSODE = Mcoef,
              diff = Mcoef-McoefGS,
              diff0 = (Mcoef==0)-(McoefGS==0),
              lambda1 = seqLambda1,
              lambda2 = seqLambda2,
              diff.range = range(Mcoef-McoefGS)))
}

