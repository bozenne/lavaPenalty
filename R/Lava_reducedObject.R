#' @title Extract variable names from a reduced latent variable model
#' @description Extract variable as in a standard lvm except that one can choose to also return the name of the linear predictors used to reduce the model
#' 
#' @param x \code{lvm}-object
#' @param lp should the name of the variables corresponding to the linear predictors be returned?
#' @param xlp should the name of the variables that the linear predictors aggregates be returned?
#' @param type slot to be return. Can be \code{"coef"}, \code{"x"}, \code{"con"}, \code{"name"}, 
#' 
#' @details lp returns all the linear predictors of the \code{lvm}-object. 
#' The other functions plays the same role as those defined in the lava package.
#' 
#' @examples 
#' 
#' ## regression
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y", reduce = TRUE)
#' vars(m)
#' vars(m, lp = FALSE)
#' vars(m, lp = FALSE, xlp = TRUE)
#' 
#' lp(m)
#' 
#' exogenous(m)
#' exogenous(m, lp = FALSE)
#' exogenous(m, xlp = TRUE)
#' 
#' endogenous(m)
#' endogenous(m, lp = TRUE) # should not change
#' endogenous(m, xlp = TRUE) # should not change
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
#' 
#' exogenous(m)
#' exogenous(m, lp = TRUE)
#' 
#' endogenous(m)
#' endogenous(m, lp = TRUE) # should not change

`lp` <- function(x,...) UseMethod("lp")

lp.lvm.reduced <- function(x, type = "name", ...){
  if(type %in% names(x$lp[[1]]) == FALSE){
    stop("type ",type," is not valid \n",
         "valid types: \"",paste(names(x$lp[[1]]), collapse = "\" \""),"\" \n")
  }
  return(unname(unlist(lapply(x$lp, function(x)x[[type]]))))
}

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

exogenous.lvm.reduced <- function(x, lp = TRUE, xlp = FALSE, ...){
  
  if(xlp){
    hiddenX <- lp(x, type = "x")
  }else{
    hiddenX <- NULL
  }
  if(lp){
    names.lp <- NULL
  }else{
    names.lp <- lp(x,...)
  }
  
  class(x) <- setdiff(class(x),"lvm.reduced")
  allExo <- unique(c(exogenous(x, ...), hiddenX))
  
  return(setdiff(allExo, names.lp))
}


# apply reduction to all equation of the LVM 
`reduce` <-
  function(x,...) UseMethod("reduce")

reduce.lvm <- function(object, response = NULL, clean.exo = TRUE){
  col.reg <- apply(object$index$Jy, 1, function(x){which(x==1)})
  cov <-  apply(object$index$A[exogenous(object),col.reg], 2,  function(x){names(which(x==1))})
  if(is.null(response)){
    response <- names(cov)#vars(object)[index.reduce]
  }
  
  if(is.null(response)){
    cat("no regression model has been found in the object \n")
    return(object)
  }
  
  for(iterR in 1:length(col.reg)){
    name.response <- response[iterR]
    name.cov <- cov[[iterR]]
    
    ## can be problematic as we don't know about "additive" or other possibly relevant arguments
    for(iter in paste(name.response,name.cov,sep ="~")){
      cancel(object) <- as.formula(iter)
    }
    
    object <- regression.lvm(object, to = name.response, from = name.cov, reduce = TRUE)
  }
  
  if(clean.exo){
    indexClean <- which(rowSums(object$index$A[object$exogenous,]!=0)==0)
    kill(object) <- object$exogenous[indexClean]
  }
  
  return(object)
}

# NOTE: the user should not be able to remove the intercept from a model that have been reduced

# add reduce argument
regression.lvm <- function(object = lvm(), to, from, fn = NA, silent = lava.options()$silent, 
                           additive = TRUE, reduce = FALSE, y, x, value, ...){
  
  if (!missing(y)) {
    if (inherits(y, "formula")) 
      y <- all.vars(y)
    to <- y
  }
  if (!missing(x)) {
    if (inherits(x, "formula")) 
      x <- all.vars(x)
    from <- x
  }
  if (!additive) {
    if (!inherits(to, "formula")) 
      to <- toformula(to, from)
    x <- attributes(getoutcome(to))$x
    K <- object$attributes$nordinal[x]
    if (is.null(K) || is.na(K)) {
      K <- list(...)$K
      if (is.null(K)) 
        stop("Supply number of categories, K (or use method 'categorical' before calling 'regression').")
      object <- categorical(object, x, ...)
    }
    dots <- list(...)
    dots$K <- K
    dots$x <- object
    dots$formula <- to
    dots$regr.only <- TRUE
    object <- do.call("categorical", dots)
    return(object)
  }
  if (missing(to)) {
    return(regfix(object))
  }
  if (inherits(to, "formula")) {
    if (!missing(value)) {
      regression(object, to, silent = silent, ...) <- value
    }
    else {
      regression(object, silent = silent, ...) <- to
    }
    object$parpos <- NULL
    return(object)
  }
  if (is.list(to)) {
    for (t in to) regression(object, silent = silent, ...) <- t
    object$parpos <- NULL
    return(object)
  }
  if(reduce){
    allCoef <- coef(regression(object, to = to, from = from, silent = silent, reduce = FALSE))
    
    if("lp" %in% names(object) == FALSE){object$lp <- list()}
    for(iterR in to){ ### I don't know where to get the constrains
      regression(object, to = iterR, from = paste0("LP",iterR), silent = silent) <- 1
      namesCoef <- paste(iterR, from, sep = "~")
      if(iterR %in% names(object$lp)){ 
        object$lp[[iterR]]$coef <- c(object$lp[[iterR]]$coef, namesCoef)
        # object$lp[[iterR]]$indexCoef <- match(object$lp[[iterR]]$indexCoef, coef(object))
        object$lp[[iterR]]$x <- c(object$lp[[iterR]]$x, from)
        object$lp[[iterR]]$con <- c(object$lp[[iterR]]$con, rep(NA,length(from)))
      }else{
        object$lp[[iterR]]$coef <- namesCoef
        # object$lp[[iterR]]$indexCoef <- match(namesCoef, coef(object))
        object$lp[[iterR]]$x <- from
        object$lp[[iterR]]$con <- rep(NA,length(from))
        object$lp[[iterR]]$name <- paste0("LP",iterR)
      }
    }
    parameter(object) <- setdiff(allCoef,coef(object))
    class(object) <- append("lvm.reduced", class(object))
    return(object)
  }
  
  sx <- strsplit(from, "@")
  xx <- sapply(sx, FUN = function(i) i[1])
  ps <- sapply(sx, FUN = function(i) i[2])
  sx <- strsplit(xx, "$", fixed = TRUE)
  xs <- sapply(sx, FUN = function(i) i[1])
  fix <- as.numeric(sapply(sx, FUN = function(i) i[2]))
  allv <- index(object)$vars
  object <- addvar(object, c(to, xs), silent = silent, reindex = FALSE)
  for (i in to) for (j in xs) {
    object$M[j, i] <- 1
    if (!is.na(fn)) 
      functional(object, j, i) <- fn
  }
  if (lava.options()$exogenous) {
    newexo <- setdiff(xs, c(to, allv))
    exo <- exogenous(object)
    if (length(newexo) > 0) 
      exo <- unique(c(exo, newexo))
    exogenous(object) <- setdiff(exo, to)
  }
  if (lava.options()$debug) {
    print(object$fix)
  }
  object$fix[xs, to] <- fix
  object$par[xs, to] <- ps
  object$parpos <- NULL
  index(object) <- reindex(object)
  
  
  return(object)
}

`cancelLP` <-
  function(x,...) UseMethod("cancelLP")

cancelLP.lvm.reduced  <- function(x, expar = TRUE){
  
  coefLP <- lp(x, type = "coef")
  
  rmvar(x) <- lp(x) 
  
  if(expar){
  x$expar <- x$expar[setdiff(names(x$expar),coefLP)]
  x$exfix <- x$exfix[setdiff(names(x$exfix),coefLP)]
  if(length(x$exfix)==0){x$exfix <- NULL}
  x$attributes$parameter <- x$attributes$parameter[setdiff(names(x$attributes$parameter),coefLP)]
  x$index$npar.ex <- x$index$npar.ex[setdiff(names(x$index$npar.ex),coefLP)]
  x$index$parval <- x$index$parval[setdiff(names(x$index$parval),coefLP)]
  }
  
  class(x) <- setdiff(class(x),"lvm.reduced")
  
  return(x)
}

