gaussianReduced_objective.lvm <- function(x, p, data, S, mu, n, implementation = "lava", ...)  {
 
  names.pReduced <- names(p)
  n.lp <- length(x$lp)
  
  for(iterLP in 1:n.lp){
    ## extract coefficients according to constrains
    b <- p[x$lp[[iterLP]]$coef] # x$lp$y1$con
    names.pReduced <- setdiff(names.pReduced, x$lp[[iterLP]]$coef)
    
    ## extract data
    form <- as.formula(paste0("~0+",paste(x$lp[[iterLP]]$x,collapse = "+")))
    X <- as.matrix(model.matrix(form, data))
    ## compute linear predictor for the reduce model
    data[,x$lp[[iterLP]]$name] <- data.frame(X %*% b)
  }
  
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
    
    l <- logLik(x$mRed,data=data,p=p[names.pReduced])
    
    # logLik(ls.save$x,data=ls.save$data,p=ls.save$p[names.pReduced])
    # ls.save <<- list(x = x,data=data,p=p[names.pReduced])
    # expect_equal(as.double(l),as.double(logLik(m,data=d,p=p[coef(m)])))
    
    
  }else if(implementation == "mets"){
    l <- mets::loglikMVN(yl = data, yu = data, status = rep(1,NROW(data)), mu = mu, S = S, thres = 0)
  }
  
  return(l)
}
gaussianReduced_logLik.lvm <- gaussianReduced_objective.lvm 

gaussianReduced_gradient.lvm <- function(x, p, data, S, mu, n, indiv = FALSE, implementation = "lava", ...)  {
  
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

    s <- score(x$mRed,data=data,p=p[names.pReduced],indiv=TRUE)
  
  }else if(implementation == "mets"){
    s <- mets::scoreMVN(yl = data, yu = data, status = rep(1,NROW(data)), 
                        mu = mu, S = S, dmu = NULL, dS = NULL)
  }  
  
  names.pHidden <- setdiff(names(p),names.pReduced)
  
  ## apply chain rule
  s.hidden <- matrix(NA, ncol = length(names.pHidden), nrow = NROW(s))
  colnames(s.hidden) <- names.pHidden
  for(iterLP in 1:length(x$lp)){ 
    name.intercept <- names(x$lp)[iterLP]
    name.lp <- x$lp[[iterLP]]$name
    
    ## extract data
    form <- as.formula(paste0("~0+",paste(x$lp[[iterLP]]$x,collapse = "+")))
    X <- as.matrix(model.matrix(form, data))
    
    dlp <- X # - what about dB/db in presence of constrains
    s.hidden[,x$lp[[iterLP]]$coef] <- apply(dlp, 2, function(j){j*s[,match(name.intercept,names.pReduced)]})
  }
  
  s.all <- cbind(s,s.hidden)[,names(p),drop = FALSE]
        
  if(indiv == FALSE){
    s.all <- apply(s.all,2,sum)
    # expect_equal(as.double(s.all[coef(m)]),as.double(score(m,data=d,p=p[coef(m)],indiv=indiv)))  
  }else{
    # expect_equal(s.all,score(m,data=d,p=p[coef(m)],indiv=indiv)[,coef(mR)])  
  }
  # derivSave <<- s.all
  return(s.all) # sum or mean
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
  return(elvm)
}

gaussianReduced_hessian.lvm <- function(x,p,...) {
  
  myg <- function(p1) gaussianReduced_gradient.lvm(x,p=p1,...)
  I <- numDeriv::jacobian(myg,p, method = "simple")
  return( (I+t(I))/2 )
  
}

### do we need a new estimate function ?????