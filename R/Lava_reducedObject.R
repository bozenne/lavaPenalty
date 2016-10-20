#' @title Extract variables from a reduced latent variable model
#' @name getReduce
#' @description Extract variables as in a standard lvm but including or not the variables that compose the linear predictor
#' 
#' @param x \code{lvm}-object
#' @param lp should the name of the linear predictors be returned?
#' @param xlp should the name of the variables that the linear predictors aggregates be returned?
#' @param type slot to be return. Can be \code{"link"}, \code{"x"}, \code{"con"}, \code{"name"}, 
#' @param index which linear predictor to consider, 
#' @param format should the results be kept as the list or returned as a single vector, 
#' 
#' @details lp returns all the linear predictors of the \code{lvm}-object. 
#' The other functions plays the same role as those defined in the lava package.
#' 
#' @examples  
#' ## regression
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' vars(m)
#' vars(m, lp = TRUE)
#' vars(m, lp = FALSE, xlp = FALSE)
#' 
#' lp(m)
#' 
#' exogenous(m)
#' exogenous(m, lp = TRUE)
#' exogenous(m, xlp = FALSE)
#' 
#' endogenous(m)
#' endogenous(m, lp = TRUE) # should not change
#' endogenous(m, xlp = FALSE) # should not change
#' 
#' coef(m)
#' 
#' ## lvm
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y1", reduce = TRUE)
#' m <- regression(m, x=paste0("x",51:150),y="y2", reduce = TRUE)
#' covariance(m) <- y1~y2
#' 
#' vars(m)
#' vars(m, lp = TRUE)
#' 
#' lp(m)
#' lp(m, type = "x", format = "list")
#' lp(m, index = 1)
#' 
#' exogenous(m)
#' exogenous(m, lp = TRUE)
#' 
#' endogenous(m)
#' endogenous(m, lp = TRUE) # should not change
#' @export
`lp` <- function(x,...) UseMethod("lp")

#' @rdname getReduce 
lp.lvm.reduced <- function(x, type = "name", lp = NULL, format = "vector", ...){
  
  if(length(x$lp)==0){return(NULL)} 
  validNames <- c("link","con","name","x") # names(x$lp[[1]])
  
  ## type
  if(is.null(type)){
    type <- validNames
    size <- FALSE
    format <- "list2"
  }else if(length(type) == 1 && type == "endogeneous"){
    return(names(x$lp))
  }else if(length(type) == 1 && type == "n.link"){
    type <- "link"
    size <- TRUE
    format <- "list"
  }else{
    if(type %in% validNames == FALSE){
      stop("type ",type," is not valid \n",
           "valid types: \"",paste(validNames, collapse = "\" \""),"\" \n")
    }
    size <- FALSE
  }
  
  ## format
  validFormat <- c("vector","list","list2")
  if(format %in% validFormat == FALSE){
    stop("format ",format," is not valid \n",
         "format must be on of : \"",paste(validFormat, collapse = " ")," \n")
  }
  if(format != "list2" && length(type)>1){
    stop("format must be \"list2\" when length(type) is not one \n",
         "length(type): ",length(type)," \n")
  }
  
  ## select lp
  if(is.null(lp)){
    lp <- seq_len(length(x$lp))
  }else if(is.numeric(lp)){
    vec <- seq_len(length(x$lp))
    if(any(lp %in% vec == FALSE)){
      stop("lp ",paste(lp, collapse = " ")," is not valid \n",
           "if numeric lp must be in: \"",paste(vec, collapse = " ")," \n")
    }
  }else if(is.character(lp)){
    vec <- unlist(lapply(x$lp, function(x)x[["name"]]))
    if(any(lp %in% vec == FALSE)){
      stop("lp ",paste(lp, collapse = " ")," is not valid \n",
           "if character lp must be in: \"",paste(vec, collapse = "\" \""),"\" \n")
    }
    lp <- match(lp, vec)
  }else{
    stop("lp must be a numeric or character vector \n")
  }
  
  ## extract 
  if(format == "list"){
    res <- lapply(x$lp[lp], function(x)x[[type]])
    
    if(size){
      res <- unlist(lapply(res, function(x) length(x)))
      names(res) <- unlist(lapply(x$lp[lp], function(x)x[["name"]]))
    }
    
  }else{
    res <- lapply(x$lp[lp], function(x)x[type])
  }
  
  ## export
  if(format == "vector"){
    res <- unlist(res)
  }
  
  return(res)
  
}

#' @rdname getReduce 
vars.lvm.reduced <- function(x, lp = TRUE, xlp = FALSE, ...){
  
  if(xlp){
    hiddenX <- lp(x, type = "x")
  }else{
    hiddenX <- NULL
  }
  if(lp){
    names.lp <- NULL
  }else{
    names.lp <- lp(x, type = "name", ...)
  }
  
  class(x) <- setdiff(class(x),"lvm.reduced")
  allVars <- unique(c(vars(x), hiddenX))
  
  return(setdiff(allVars,names.lp))
}

#' @rdname getReduce 
exogenous.lvm.reduced <- function(x, lp = TRUE, xlp = FALSE, ...){
  
  if(xlp){
    hiddenX <- lp(x, type = "x")
  }else{
    hiddenX <- NULL
  }
  if(lp){
    names.lp <- NULL
  }else{
    names.lp <- lp(x, type = "name")
  }
  
  class(x) <- setdiff(class(x),"lvm.reduced")
  allExo <- unique(c(exogenous(x), hiddenX))
  
  return(setdiff(allExo, names.lp))
}

#' @rdname getReduce 
endogenous.lvm.reduced <- function(x, top = FALSE, latent = FALSE, ...){
  
  observed <- manifest(x, lp = FALSE, xlp = FALSE)
  if (latent) observed <- vars(x, lp = FALSE, xlp = FALSE)
  if (top) {
    M <- x$M
    res <- c()
    for (i in observed)
      if (!any(M[i,]==1))
        res <- c(res, i)
    return(res)
  }
  exo <- exogenous(x, lp = FALSE, xlp = FALSE)
  return(setdiff(observed,exo))
  
}

#' @rdname getReduce 
manifest.lvm.reduced <- function(x, lp = TRUE, xlp = FALSE, ...) {
  
  vars <- vars(x, lp = lp, xlp = xlp)
  if (length(vars)>0)
    setdiff(vars,latent(x))
  else
    NULL
}

#' @title Reduce latent variable model
#' @name reduce
#' @description Aggregate exogeneous variables into a linear predictor
#' 
#' @param x \code{lvm}-object
#' @param endo the endogeneous variables for which the related exogeneous variables should be aggregated
#' @param rm.exo should the exogeneous variables be remove from the object
#' 
#' @examples 
#' ## regression
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' 
#' m
#' reduce(m)
#' coef(reduce(m))
#' reduce(m, rm.exo = FALSE)
#' 
#' ## lvm
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y1")
#' m <- regression(m, x=paste0("x",51:150),y="y2")
#' covariance(m) <- y1~y2
#' 
#' m
#' reduce(m)
#' coef(reduce(m))
#' reduce(m, rm.exo = FALSE)
#' 
#' ## with penalty
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m)
#' reduce(pm)
#' reduce(pm, link = y~x1+x2+x3)
#' 
`reduce` <-
  function(x,...) UseMethod("reduce")

#' @rdname getReduce 
reduce.lvm <- function(object, link = NULL, endo = NULL, rm.exo = TRUE){
  
  if(!is.null(link)){ # reduce specific links
    
    if("formula" %in% class(link)){
      ls.link <- initVar_link(link, repVar1 = TRUE)
      vec.endo <- ls.link$var1
      vec.exo <- ls.link$var2
    }else{
      ls.link <- lapply(link, initVar_link, repVar1 = TRUE)
      vec.endo <- unlist(lapply(ls.link, "[[", 1))
      vec.exo <- unlist(lapply(ls.link, "[[", 2))
    }
    
    exo <- tapply(vec.exo, vec.endo, list)
    endo <- names(exo)
    
  }else{ # reduce the linear predictor of specific endogeneous variables
    col.reg <- apply(object$index$Jy, 1, function(x){which(x==1)}) # find endo having exo
    
    if(length(col.reg)==1){ # find exo corresponding to endo
      exo <- list(names(which(object$index$A[exogenous(object),col.reg]==1)))
      if(length(exo[[1]])>0){names(exo) <- colnames(object$index$A)[col.reg]}
    }else{
      exo <-  apply(object$index$A[exogenous(object),col.reg, drop = FALSE], 2,  function(x){names(which(x==1))})
    }
    
    if(is.null(endo)){
      if(is.null(names(exo))){
        cat("no regression model has been found in the object \n")
        return(object)
      }else{
        endo <- names(exo)#vars(object)[index.reduce]
      }
    }
    
  }
  
  n.endo <- length(endo)
  
  for(iterR in 1:n.endo){
    name.endo <- endo[iterR]
    name.exo <- exo[[iterR]]
    
    ## can be problematic as we don't know about "additive" or other possibly relevant arguments
    f <- as.formula(paste(name.endo,"~",paste(name.exo, collapse = "+")))
    cancel(object) <- f
    
    object <- regression.lvm(object, to = name.endo, from = name.exo, reduce = TRUE)
  }
  
  if(rm.exo){
    indexClean <- which(rowSums(object$index$A[object$exogenous,]!=0)==0)
    kill(object) <- object$exogenous[indexClean]
  }
  
  return(object)
}

#' @rdname getReduce 
reduce.plvm <- function(object, link = NULL, rm.exo = TRUE, ...){
  
  if(is.null(link)){
    link <- penalty(object, type = "link", nuclear = FALSE) 
  }
  object <- reduce.lvm(object, link = link, rm.exo = rm.exo, ...)
  
  return(object)
}

#' @title Update of the linear predictor
#' @description Update one or several linear predictors
#' 
#' @param x \code{lvm}-object
#' @param lp the name of the linear predictors to be updated
#' @param value the value that will be allocated
#' 
#' @examples 
#' 
#' ## regresssion
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' vars(m)
#' 
#' newLP <- lp(m, type = NULL)[[1]]
#' newLP$link <- newLP$link[1:3]
#' newLP$con <- newLP$con[1:3]
#' newLP$x <- newLP$x[1:3]
#' 
#' lp(m, lp = 1) <- newLP
#' lp(m, type = NULL)
#' 
`lp<-` <- function(x,...) UseMethod("lp<-")

`lp<-.lvm.reduced` <- function(x, lp = NULL, value){
  
  if(is.null(lp)){
    
    x$lp <- value
    
  }else{
    
    ## valid LP
    if(is.numeric(lp)){
      vec <- seq_len(length(x$lp))
      if(any(lp %in% vec == FALSE)){
        stop("lp ",paste(lp, collapse = " ")," is not valid \n",
             "if numeric lp must be in: \"",paste(vec, collapse = " ")," \n")
      }
    }else if(is.character(lp)){
      vec <- unlist(lapply(x$lp, function(x)x[["name"]]))
      if(any(lp %in% vec == FALSE)){
        stop("lp ",paste(lp, collapse = " ")," is not valid \n",
             "if character lp must be in: \"",paste(vec, collapse = " ")," \n")
      }
      lp <- match(lp, vec)
    }
    
    if(length(lp) == 1 && length(value) == 4 && all(names(value) == c("link","con","name","x"))){
      value <- list(value)
    }
    
    x$lp[lp] <- value
    
  }
  
  return(x)
  
}

#' @title Remove variable from the linear predictor
#' @description Remove one or several variable from the linear predictor
#' @name cancelLP
#' 
#' @param x \code{lvm}-object
#' @param lp should the name of the variables corresponding to the linear predictors be returned?
#' @param expar should the external parameters be also removed
#' 
#' @examples 
#' 
#' ## regresssion
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' vars(m)
#' 
#' cancelLP(m)
#' cancelLP(m, lp = "LPy", link = c("y~x1","y~x2"))
#' 
#' ## with penalization
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m, reduce = TRUE)
#' cancelLP(pm)
#' coef(cancelLP(pm))


#' @export
#' @rdname cancelLP
`cancelLP` <-
  function(x,...) UseMethod("cancelLP")

#' @export
#' @rdname cancelLP
`cancelLP<-` <- function (x, ..., value) {
  UseMethod("cancelLP<-", x)
}

`cancelLP<-.lvm.reduced` <- function (x, ..., value) {
  return(cancelLP(x, link = value, ...))
}

cancelLP.lvm.reduced  <- function(x, link = NULL, lp = NULL, expar = TRUE, simplify = TRUE){
  
  if(is.null(lp)){
    lp <- lp(x) 
  }else if(!is.list(lp)){
    lp <- as.list(lp)
  }
  n.lp <- length(lp)
  names.lp <- unlist(lp)
  
  if(is.null(link)){
    link <- lp(x, type = "link", lp = unlist(lp), format = "list")
  }else if(!is.list(link)){
    link <- lapply(1:n.lp, function(x) link)
    names(link) <- names.lp
  }
  
  ## remove external parameters from the LVM 
  if(expar){
    parameter(x, remove = TRUE) <- unlist(link)
    # x$expar <- x$expar[setdiff(names(x$expar),coefLP)]
    # x$exfix <- x$exfix[setdiff(names(x$exfix),coefLP)]
    # if(length(x$exfix)==0){x$exfix <- NULL}
    # x$attributes$parameter <- x$attributes$parameter[setdiff(names(x$attributes$parameter),coefLP)]
    # x$index$npar.ex <- x$index$npar.ex[setdiff(names(x$index$npar.ex),coefLP)]
    # x$index$parval <- x$index$parval[setdiff(names(x$index$parval),coefLP)]
  }
  
  ## removing variable in the LP
  for(iterLP in n.lp){
    newlp <- lp(x, type = NULL, lp = names.lp[iterLP])[[1]]
    
    indexRM <- which(newlp$link %in% link[[iterLP]])
    
    newlp$link <- newlp$link[-indexRM]
    newlp$con <- newlp$con[-indexRM]
    newlp$x <- newlp$x[-indexRM]
    
    lp(x, lp = names.lp[iterLP]) <- newlp
  }
  
  ## remove penalization
  if("plvm" %in% class(x)){
    x <- cancelPenalty(x, link = unlist(link))
  }
  
  ## remove empty LP
  n.link <- lp(x, type = "n.link")
  if(any(n.link==0)){
    rmvar(x) <- names(n.link)[which(n.link==0)]
  }
  
  ## update class
  if(simplify && length(lp(x, type = "link", format = "vector")) == 0 ){
    class(x) <- setdiff(class(x), "lvm.reduced")
  }
  
  return(x)
}

