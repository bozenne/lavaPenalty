gaussianLP_method.lvm <- "nlminb2"

gaussianLP_objective.lvm <- function(x, p, data, ...){ 
  l <- gaussianLP_logLik.lvm(x, p=p, data=data, ...)
  return(-l)
}


gaussianLP_logLik.lvm <- function(object, p, data, ...)  {
  
  dataLP <- calcLP.lvm(object, p = p, data = data, 
                       lp.x = lp(object, type = "x", format = "list"), 
                       lp.link = lp(object, type = "link", format = "list"), 
                       lp.endo = lp(object, type = "endogeneous"))
  
  ## from normal_objective.lvm
  y <- lava::index(object)$endogenous
  ord <- lava::ordinal(object)
  status <- rep(0,length(y))
  bin <- tryCatch(match(do.call("binary",list(x=object)),y),error=function(x) NULL)
  status[match(ord,y)] <- 2
  
  mu <- predict(object,data=dataLP,p=p) ## can be optimized
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

gaussianLP_gradient.lvm <- function(x, p, data, ...){
  val <- -gaussianLP_score.lvm(x, p = p, data = data, ...)
  if (!is.null(nrow(val))) {
    val <- colSums(val)
  }
  # val <- unname(val)
  return(val)
}

#' @title Compute the score for a reduced lvm model
#'
#' @details this function assumes that the external parameters in p are at the end of the vector
#' 
gaussianLP_score.lvm <- function(x, p, data, indiv = FALSE, ...)  {
  
  if(is.matrix(p)){
    p <- as.double(p)
  }
  lp.endo <- lp(x, type = "endogeneous")
  lp.link <- lp(x, type = "link", format = "list")
  lp.x <- lp(x, type = "x", format = "list")
  dataLP <- calcLP.lvm(x, p = p, data = data, 
                       lp.x = lp.x, lp.link = lp.link, lp.endo = lp.endo)
  
  ## from normal_gradient.lvm
  M <- moments(x,p)
  Y <- as.matrix(dataLP[,manifest(x)])
  D <- deriv.lvm(x,p=p)
  mu <- M$xi%x%rep(1,nrow(Y))
  
  s <- mets::scoreMVN(Y,mu,M$C,D$dxi,D$dS)
  colnames(s) <- coef(x)#c(name.intercept,name.regression,name.covariance)
  
  ## apply chain rule
  n.lp <- length(x$lp)
  
  for(iterLP in 1:n.lp){ 
    name.intercept <- lp.endo[iterLP]
    
    ## extract data
    X <-  data[,lp.x[[iterLP]],drop = FALSE]
    dlp <- X # - what about dB/db in presence of constrains
    s[,lp.link[[iterLP]]] <- apply(dlp, 2, function(j){j*s[,name.intercept]})
  }
  
  if(indiv == FALSE){
    s <- rbind(apply(s,2,sum))
  }
  print(s)
  return(s) 
}

gaussianLP_hessian.lvm <- function(x,p,n,type,...) {
  dots <- list(...);
  
  if(type == "E"){
    S <- -gaussianLP_score.lvm(x,p=p,data=dots$data,indiv = TRUE)
    I <- t(S)%*%S
    attributes(I)$grad <- colSums(S)
    return(I)
  }else if(type=="num"){
    myg <- function(p1){ -gaussianLP_score.lvm(x,p=p1,n=n,data=dots$data,indiv=FALSE) }
    I <- numDeriv::jacobian(myg,p) # , method = "simple"
    I <- (I+t(I))/2
    return( I )
  }else if(type == "information"){ ## true part
    
    lp.endo <- lp(x, type = "endogeneous")
    lp.link <- lp(x, type = "link", format = "list")
    lp.x <- lp(x, type = "x", format = "list")
    dataLP <- calcLP.lvm(x, p = p, data = data, 
                         lp.x = lp.x, lp.link = lp.link, lp.endo = lp.endo)
    
    ## direct
    I <- information(x=x, p=p, n=n, data = dataLP)
    return(I)
    
    ## indirect
    
    names.pHidden <- setdiff(names(p),resLP$names$reduced)
    
    I.hiddenW <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(resLP$names$reduced))
    I.hiddenB <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(names.pHidden)) 
    
    I.all <- rbind(cbind(I, I.hiddenW), cbind(t(I.hiddenW),I.hiddenB))
    
    return(I.all)
    stop("gaussianLP_hessian.lvm not ready \n")
    # ## apply chain rule
    # I.hiddenW <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(s))
    # I.hiddenB <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(s)) 
    
    
  }
  
}

gaussian1LP_logLik.lvm <- gaussianLP_logLik.lvm
gaussian1LP_objective.lvm <- gaussianLP_objective.lvm
gaussian1LP_score.lvm <- gaussianLP_score.lvm
gaussian1LP_gradient.lvm <- gaussianLP_gradient.lvm
gaussian1LP_hessian.lvm <- function(x, type, ...){
  gaussianLP_hessian.lvm(x, type = "num", ...)
}

gaussian2LP_logLik.lvm <- gaussianLP_logLik.lvm
gaussian2LP_objective.lvm <- gaussianLP_objective.lvm
gaussian2LP_score.lvm <- gaussianLP_score.lvm
gaussian2LP_gradient.lvm <- gaussianLP_gradient.lvm
gaussian2LP_hessian.lvm <- function(x, type, ...){
  gaussianLP_hessian.lvm(x, type = "E", ...)
}



#' @title Compute the linear predictor
#' @description Compute the value of the linear predictors of a LVM and store it into the dataset
calcLP.lvm <- function(x, p, data,
                       lp.x, lp.link, lp.endo){
  
  pext <- modelVar(x,p)$e # redefine p with names (only external parameters)
  lp.name <- lp(x, type = "name")
  n.lp <- length(lp.name)
  
  for(iterLP in 1:n.lp){
    ## extract coefficients according to constrains
    b <- pext[lp.link[[iterLP]]] # x$lp$y1$con
    
    ## compute linear predictor for the reduce model
    # form <- as.formula(paste0("~0+",paste(lp.x,collapse = "+")))
    # as.matrix(model.matrix(form, data))
    # if(is.data.table(data)){
    #   data[,lp.name[[iterLP]] := as.matrix(.SD) %*% b, with = FALSE, .SDcols = lp.x[[iterLP]]] 
    # }else{
      X <- as.matrix(data[,lp.x[[iterLP]],drop = FALSE]) 
      data[,lp.name[iterLP]] <- X %*% b
    # }
   
  }
  return(data)
}
