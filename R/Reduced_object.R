#' @title Reduce latent variable model
#' @name reduce
#' @description Aggregate exogeneous variables into a linear predictor
#' 
#' @param x \code{lvm}-object
#' @param endo the endogeneous variables for which the related exogeneous variables should be aggregated
#' @param rm.exo should the exogeneous variables be remove from the object
#' 
#' @examples  
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m)
#' reduce(pm)
#' reduce(pm, link = y~x1+x2+x3)
#'

#' @rdname getReduce 
#' @export
reduce.plvm <- function(object, link = NULL, rm.exo = TRUE, ...){
  
  if(is.null(link)){
    link <- penalty(object, type = "link", nuclear = FALSE) 
  }
  object <- reduce.lvm(object, link = link, rm.exo = rm.exo, ...)
  
  return(object)
}



#' @title Remove variable from the linear predictor
#' @description Remove one or several variable from the linear predictor
#' @name cancelLP
#' 
#' @param x \code{lvm}-object
#' @param link,value the name of the links that should be removed
#' @param lp the name of the linear predictors that should be removed
#' @param expar should the external parameters be also removed
#' @param restaure should the link be kept while removing the linear predictor
#' @param simplify should the class \code{lvm.reduced} be removed is from the \code{lvm}-object if it contains no LP.
#' @param ... argument passed to \code{cancelLP.lvm.reduced}
#' 
#' @examples 
#' 
#' @export
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
  for(iterLP in 1:n.lp){
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

