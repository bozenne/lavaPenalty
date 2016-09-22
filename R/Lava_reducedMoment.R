gaussianReduced_objective.lvm <- function(x, p, data, S, mu, n, ...){ # NOT CORRECT

  k <- length(index(x)$endogenous)#length(index(x)$manifest)-length(x$lp)
  correction <- (n*k)/2*log(2*base::pi)
  
  ll <- as.double(gaussianReduced_logLik.lvm(x, p=p, data=data, S=S, mu=mu, n=n,...))
  ll <- ll + correction
  return(-ll)
}

gaussianReduced_logLik.lvm <- function(x, p, data, S, mu, n, implementation = "lava", ...)  {
 
  resLP <- calcLP.lvm(x, p = p, data = data)
  
  ## compute the log-lik of the reduced model
  if(implementation == "lava"){
    
    # l <- logLik(x,data=data,p=p[names.pReduced]) [ERROR due to the additonnal parameters]
    
    # ## need to remove additional parameters
    # additionalParam <- unlist(lapply(x$lp,function(x){x$coef}))
    # index.rm <- which(x$index$eparname.all %in% additionalParam)
    # x$index$e0 <- x$index$e0[-index.rm]
    # x$index$e1 <- x$index$e0[-index.rm]
    # x$index$npar.ex <- x$index$npar.ex - length(index.rm)
    # x$index$eparname.all <- x$index$eparname.all[-index.rm]
    # x$index$eparname <- x$index$eparname[-index.rm]
    # x$index$eparname.all.idx <- x$index$eparname.all.idx[-index.rm]
    # x$index$epar <- x$index$epar[-index.rm,]
    
    l <- logLik(x$mRed,data=resLP$data,p=p[resLP$names$reduced])
    
    # logLik(ls.save$x,data=ls.save$data,p=ls.save$p[names.pReduced])
    # ls.save <<- list(x = x,data=data,p=p[names.pReduced])
    # expect_equal(as.double(l),as.double(logLik(m,data=d,p=p[coef(m)])))
    
    
  }else if(implementation == "mets"){
    l <- mets::loglikMVN(yl = resLP$data, yu = resLP$data, status = rep(1,NROW(resLP$data)), mu = mu, S = S, thres = 0)
  }
  
  return(l)
}

gaussianReduced_gradient.lvm <- function(x, p, data, S, mu, n, ...){
  
  val <- -gaussianReduced_score.lvm(x, p = p, S = S, mu = mu, n = n, 
                                    data = data, reindex = FALSE, ...)
  if (!is.null(nrow(val))) {
    val <- colSums(val)
  }
  
  return(val)
}

gaussianReduced_score.lvm <- function(x, p, data, S, mu, n, indiv = FALSE, implementation = "lava", ...)  {
  
  resLP <- calcLP.lvm(x, p = p, data = data)
  
  ## compute the log-lik of the reduced model using mets
  if(implementation == "lava"){
    
    # s <- score(x,data=data,p=p[names.pReduced],indiv=TRUE) # [ERROR due to the additonnal parameters]
    
    # ## need to remove additional variables
    # rmvar(x) <- names.varLP
    # 
    # ## need to remove additional parameters
    # additionalParam <- unlist(lapply(x$lp,function(x){x$coef}))
    # index.rm <- which(x$index$eparname.all %in% additionalParam)
    # x$index$e0 <- x$index$e0[-index.rm]
    # x$index$e1 <- x$index$e0[-index.rm]
    # x$index$npar.ex <- x$index$npar.ex - length(index.rm)
    # x$index$eparname.all <- x$index$eparname.all[-index.rm]
    # x$index$eparname <- x$index$eparname[-index.rm]
    # x$index$eparname.all.idx <- x$index$eparname.all.idx[-index.rm]
    # x$index$epar <- x$index$epar[-index.rm,]

    s <- score(x$mRed,data=resLP$data,p=p[resLP$names$reduced],indiv=TRUE)
  
  }else if(implementation == "mets"){
    s <- mets::scoreMVN(yl = resLP$data, yu = resLP$data, status = rep(1,NROW(resLP$data)), 
                        mu = mu, S = S, dmu = NULL, dS = NULL)
  }  
  
  names.pHidden <- setdiff(names(p),resLP$names$reduced)
  
  ## apply chain rule
  s.hidden <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(s))
  colnames(s.hidden) <- names.pHidden
  for(iterLP in 1:length(x$lp)){ 
    name.intercept <- names(x$lp)[iterLP]
    name.lp <- x$lp[[iterLP]]$name
    
    ## extract data
    form <- as.formula(paste0("~0+",paste(x$lp[[iterLP]]$x,collapse = "+")))
    X <- as.matrix(model.matrix(form, resLP$data))
    
    dlp <- X # - what about dB/db in presence of constrains
    s.hidden[,x$lp[[iterLP]]$coef] <- apply(dlp, 2, function(j){j*s[,match(name.intercept,resLP$names$reduced)]})
  }
  
  s.all <- cbind(s,s.hidden)[,names(p),drop = FALSE]
        
  if(indiv == FALSE){
    s.all <- apply(s.all,2,sum)
  }
  
  return(rbind(s.all))
}

gaussianReduced_hessian.lvm <- function(x,p,n,data,type,implementation = "lava",...) {
 
   if(type == "E"){
      S <- -gaussianReduced_score.lvm(x,p=p,n=n,data=data,...)
      I <- t(S)%*%S
      attributes(I)$grad <- colSums(S)
      return(I)
    }else if(type=="num"){
      myg <- function(p1) gaussianReduced_gradient.lvm(x,p=p1,n=n,data=data,indiv=FALSE,...)
      I <- numDeriv::jacobian(myg,p, method = "simple")
      I <- (I+t(I))/2
      #attributes(I)$grad <- myg(p)
      return( I )
    }else if(type == "information"){ ## true part
      
    resLP <- calcLP.lvm(x, p = p, data = data)
    
    if(implementation == "lava"){
      
      dots <- list(...)
      dots$weight <- NULL
      I <- information(x$mRed,data=resLP$data,p=p[resLP$names$reduced],indiv=TRUE)
      # do.call("information", c(list(x = x$mRed, p = p[resLP$names$reduced], n = n, data = resLP$data), dots))
      
    }else if(implementation == "mets"){
      stop("informationMVN not yet implemented in mets \n")
      
    } 
    
    names.pHidden <- setdiff(names(p),resLP$names$reduced)
    
    I.hiddenW <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(resLP$names$reduced))
    I.hiddenB <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(names.pHidden)) 
    
    I.all <- rbind(cbind(I, I.hiddenW), cbind(t(I.hiddenW),I.hiddenB))
   
    return(I.all)
    stop("gaussianReduced_hessian.lvm not ready \n")
    # ## apply chain rule
    # I.hiddenW <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(s))
    # I.hiddenB <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(s)) 
    # s.hidden <- 
    #   colnames(s.hidden) <- names.pHidden
    # for(iterLP in 1:length(x$lp)){ 
    #   name.intercept <- names(x$lp)[iterLP]
    #   name.lp <- x$lp[[iterLP]]$name
    #   
    #   ## extract data
    #   form <- as.formula(paste0("~0+",paste(x$lp[[iterLP]]$x,collapse = "+")))
    #   X <- as.matrix(model.matrix(form, resLP$data))
    #   
    #   dlp <- X # - what about dB/db in presence of constrains
    #   s.hidden[,x$lp[[iterLP]]$coef] <- apply(dlp, 2, function(j){j*s[,match(name.intercept,names.pReduced)]})
    # }
    # 
    # s.all <- cbind(s,s.hidden)[,names(p),drop = FALSE]
    # 
    # if(indiv == FALSE){
    #   s.all <- apply(s.all,2,sum)
    # }
    
  }
  
}

estimate.lvm.reduced <- function(x, data, control = list(), ...){

  ## add observations for LP
  names.LP <- unlist(lapply(x$lp, function(j){j$name}))
  data <- cbind(data,
                data.frame(matrix(0,nrow = NROW(data), ncol = length(names.LP), dimnames = list(NULL,names.LP)))
  )
  
  ## intialisation of the LP paramters
  if("start" %in% names(control) == FALSE){
    ls.coef <- lapply(names(x$lp), function(j){
      coef <- lm.fit(y = data[[j]], x = as.matrix(data[x$lp[[j]]$x]))$coefficients
      names(coef) <- paste(j,names(coef),sep="~")
      return(coef)
    })
    control$start <- unlist(ls.coef)[na.omit(match(coef(x),names(unlist(ls.coef))))]
  }
  
  elvm <- lava:::estimate.lvm(x, data = data, control = control, ...)
  # elvm <- estimate.lvm(x, data = data, control = control, ...)
  return(elvm)
}


#' @description Compute the value of the linear predictors of a LVM and store it into the dataset
calcLP.lvm <- function(x, p, data){
  
  names.pReduced <- names(p)
  names.varLP <- NULL
  n.lp <- length(x$lp)
  
  for(iterLP in 1:n.lp){
    ## extract coefficients according to constrains
    b <- p[x$lp[[iterLP]]$coef] # x$lp$y1$con
    names.pReduced <- setdiff(names.pReduced, x$lp[[iterLP]]$coef)
    
    ## extract data
    form <- as.formula(paste0("~0+",paste(x$lp[[iterLP]]$x,collapse = "+")))
    X <- as.matrix(model.matrix(form, data))
    names.varLP <- c(names.varLP, x$lp[[iterLP]]$x)
    
    ## compute linear predictor for the reduce model
    data[,x$lp[[iterLP]]$name] <- data.frame(X %*% b)
  }
  
  return(list(data = data,
              names = list(varLP = names.varLP, reduced = names.pReduced)))
}
