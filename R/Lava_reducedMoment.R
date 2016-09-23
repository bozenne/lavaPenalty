gaussianReduced_objective.lvm <- function(x, p, data, ...){ 
  
  l <- gaussianReduced_logLik.lvm(x, p=p, data=data, ...)
  return(-l)
}

gaussianReduced_logLik.lvm <- function(x, p, data, ...)  {
  
  resLP <- calcLP.lvm(x, p = p, data = data)
  
  ## from normal_objective.lvm
  y <- lava::index(x)$endogenous
  ord <- lava::ordinal(x)
  status <- rep(0,length(y))
  bin <- tryCatch(match(do.call("binary",list(x=x)),y),error=function(x) NULL)
  status[match(ord,y)] <- 2
  
  mu <- predict(x$mRed,data=resLP$data,p=p)## --mRed
  S <- attributes(mu)$cond.var
  class(mu) <- "matrix"
  thres <- matrix(0,nrow=length(y),max(1,attributes(ord)$K-1)); rownames(thres) <- y
  for (i in seq_len(length(attributes(ord)$fix))) {
    nn <- names(attributes(ord)$idx)[i]
    ii <- attributes(ord)$idx[[nn]]
    val <- (attributes(mu)$e[ii])
    thres[nn,seq_len(length(val))] <-
      cumsum(c(val[1],exp(val[-1])))
  }
  
  yl <- yu <- as.matrix(data[,y,drop=FALSE])
  
  if (!inherits(yl[1,1],c("numeric","integer","logical")) ||
      !inherits(yu[1,1],c("numeric","integer","logical")))
    stop("Unexpected data (normal_objective)")
  
  l <- mets::loglikMVN(yl = yl, yu = yu, status = status, mu = mu, S = S, thres = thres)
  l <- sum(l)
  
  return(l)
}

gaussianReduced_gradient.lvm <- function(x, p, data, ...){
  
  val <- -gaussianReduced_score.lvm(x, p = p, data = data, ...)
  if (!is.null(nrow(val))) {
    val <- colSums(val)
  }
  
  return(val)
}

gaussianReduced_score.lvm <- function(x, p, data, indiv = FALSE, ...)  {
  
  resLP <- calcLP.lvm(x, p = p, data = data)
  
  ## from normal_gradient.lvm
  D <- lava:::deriv.lvm(x$mRed,p=p) ## --mRed
  M <- moments(x$mRed,p)
  Y <- as.matrix(resLP$data[,manifest(x$mRed)])
  mu <- M$xi%x%rep(1,nrow(Y))
  s <-mets::scoreMVN(Y,mu,M$C,D$dxi,D$dS)
  colnames(s) <-  coef(x$mRed)
  
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

gaussianReduced_hessian.lvm <- function(x,p,n,type,...) {
  
  dots <- list(...);
  # dots$weight <- NULL
  # return(do.call("information", c(list(x=x,p=p,n=n),dots)))
  
  if(type == "E"){
    S <- -gaussianReduced_score.lvm(x,p=p,data=dots$data,indiv=FALSE)
    I <- t(S)%*%S
    attributes(I)$grad <- colSums(S)
    return(I)
  }else if(type=="num"){
    myg <- function(p1) gaussianReduced_gradient.lvm(x,p=p1,n=n,data=dots$data,indiv=FALSE)
    I <- numDeriv::jacobian(myg,p, method = "simple")
    I <- (I+t(I))/2

    return( I )
  }else if(type == "information"){ ## true part
    
    resLP <- calcLP.lvm(x, p = p, data = dots$data)
    
    ## direct
    I <- information(x=x, p=p, n=n, data = resLP$data)
    return(I)
    
    ## indirect
    I <- information(x$mRed,data=resLP$data,p=p[resLP$names$reduced],indiv=TRUE)
    
    names.pHidden <- setdiff(names(p),resLP$names$reduced)
    
    I.hiddenW <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(resLP$names$reduced))
    I.hiddenB <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(names.pHidden)) 
    
    I.all <- rbind(cbind(I, I.hiddenW), cbind(t(I.hiddenW),I.hiddenB))
    
    return(I.all)
    stop("gaussianReduced_hessian.lvm not ready \n")
    # ## apply chain rule
    # I.hiddenW <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(s))
    # I.hiddenB <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(s)) 
    
    
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
