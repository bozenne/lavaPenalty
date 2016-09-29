gaussianReduced_objective.lvm <- function(x, p, data, ...){ 
  l <- gaussianReduced_logLik.lvm(x, p=p, data=data, ...)
  return(-l)
}


gaussianReduced_logLik.lvm <- function(x, p, data, ...)  {
  dataLP <- calcLP.lvm(x, p = p, data = data)
  
  ## from normal_objective.lvm
  y <- lava::index(x)$endogenous
  ord <- lava::ordinal(x)
  status <- rep(0,length(y))
  bin <- tryCatch(match(do.call("binary",list(x=x)),y),error=function(x) NULL)
  status[match(ord,y)] <- 2
  
  mu <- predict(x,data=dataLP,p=p) ## can be optimized
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
  # val <- unname(val)
  return(val)
}

#' @details this function assumes that the external parameters in p are at the end of the vector
#' 
gaussianReduced_score.lvm <- function(x, p, data, indiv = FALSE, ...)  {
  
  if(is.matrix(p)){
    p <- as.double(p)
  }
  dataLP <- calcLP.lvm(x, p = p, data = data)
  
  ## from normal_gradient.lvm
  M <- moments(x,p)
  Y <- as.matrix(dataLP[,manifest(x)])
  D <- deriv.lvm(x,p=p)
  mu <- M$xi%x%rep(1,nrow(Y))
  
  
  s <- mets::scoreMVN(Y,mu,M$C,D$dxi,D$dS)
  # print(head(s))
  # name.intercept <- intersect(names(M$v),coef(x))
  # name.regression <- names(M$e)
  # name.covariance <- setdiff(coef(x), c(name.intercept,name.regression))
  colnames(s) <- coef(x)#c(name.intercept,name.regression,name.covariance)
  
  ## apply chain rule
  n.lp <- length(x$lp)
  for(iterLP in 1:n.lp){ 
    name.intercept <- names(x$lp)[iterLP]
    name.lp <- x$lp[[iterLP]]$name
    
    ## extract data
    form <- as.formula(paste0("~0+",paste(x$lp[[iterLP]]$x,collapse = "+")))
    X <- as.matrix(model.matrix(form, dataLP))
    
    dlp <- X # - what about dB/db in presence of constrains
    s[,x$lp[[iterLP]]$coef] <- apply(dlp, 2, function(j){j*s[,name.intercept]})
  }
  
  if(indiv == FALSE){
    s <- rbind(apply(s,2,sum))
  }
  
  return(s) 
}

gaussianReduced_hessian.lvm <- function(x,p,n,type,...) {
  dots <- list(...);

  if(type == "E"){
    S <- -gaussianReduced_score.lvm(x,p=p,data=dots$data,indiv = TRUE)
    I <- t(S)%*%S
    attributes(I)$grad <- colSums(S)
    return(I)
  }else if(type=="num"){
    myg <- function(p1){ -gaussianReduced_score.lvm(x,p=p1,n=n,data=dots$data,indiv=FALSE) }
    I <- numDeriv::jacobian(myg,p) # , method = "simple"
    I <- (I+t(I))/2
    return( I )
  }else if(type == "information"){ ## true part
    
    dataLP <- calcLP.lvm(x, p = p, data = dots$data)
    
    ## direct
    I <- information(x=x, p=p, n=n, data = dataLP)
    return(I)
    
    ## indirect
    
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

gaussianReduced1_logLik.lvm <- gaussianReduced_logLik.lvm
gaussianReduced1_objective.lvm <- gaussianReduced_objective.lvm
gaussianReduced1_score.lvm <- gaussianReduced_score.lvm
gaussianReduced1_gradient.lvm <- gaussianReduced_gradient.lvm
gaussianReduced1_hessian.lvm <- function(x, type, ...){
  gaussianReduced_hessian.lvm(x, type = "num", ...)
}

gaussianReduced2_logLik.lvm <- gaussianReduced_logLik.lvm
gaussianReduced2_objective.lvm <- gaussianReduced_objective.lvm
gaussianReduced2_score.lvm <- gaussianReduced_score.lvm
gaussianReduced2_gradient.lvm <- gaussianReduced_gradient.lvm
gaussianReduced2_hessian.lvm <- function(x, type, ...){
  gaussianReduced_hessian.lvm(x, type = "E", ...)
}


estimate.lvm.reduced <- function(x, data, control = list(), estimator = "gaussianReduced", ...){
  
  ## add observations for LP
  names.LP <- unlist(lapply(x$lp, function(j){j$name}))
  data <- cbind(data,
                data.frame(matrix(0,nrow = NROW(data), ncol = length(names.LP), dimnames = list(NULL,names.LP)))
  )
  
  ## intialisation of the LP paramters
  if("start" %in% names(control) == FALSE){
    
    x0 <- cancelLP(x)
    startLVM <- coef(lava:::estimate.lvm(x0, data))
    
    ls.coef <- lapply(names(x$lp), function(j){
      coef <- lm.fit(y = data[[j]], x = as.matrix(data[x$lp[[j]]$x]))$coefficients
      names(coef) <- paste(j,names(coef),sep="~")
      return(coef)
    })
    startLP <- unlist(ls.coef)[na.omit(match(coef(x),names(unlist(ls.coef))))]
    
    
    control$start <- setNames(rep(NA, length = length(coef(x))), coef(x))
    control$start[names(startLP)] <- startLP
    control$start[names(startLVM)] <- startLVM
    control$start <- control$start[which(!is.na(control$start))]
  }
  
  print(control$start)
  
  # elvm <- lava:::estimate.lvm(x, data = data, control = control, ...)
  elvm <- estimate.lvm(x, data = data, control = control, estimator = estimator, ...)
  return(elvm)
}

#' @description Compute the value of the linear predictors of a LVM and store it into the dataset
calcLP.lvm <- function(x, p, data){
  
  pext <- modelVar(x,p)$e # redefine p with names (only external parameters)
  n.lp <- length(x$lp)
  for(iterLP in 1:n.lp){
    ## extract coefficients according to constrains
    b <- pext[x$lp[[iterLP]]$coef] # x$lp$y1$con
   
    ## extract data
    form <- as.formula(paste0("~0+",paste(x$lp[[iterLP]]$x,collapse = "+")))
    X <- as.matrix(model.matrix(form, data))
    
    ## compute linear predictor for the reduce model
    data[,x$lp[[iterLP]]$name] <- data.frame(X %*% b)
  }
  
  
  return(data)
}
