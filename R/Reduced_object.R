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
  object <- lava.reduce:::reduce.lvm(object, link = link, rm.exo = rm.exo, ...)
  
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
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m)
#' rpm <- reduce(pm)
#' 
#' cancelLP(rpm)
#' 
#' reduce(pm, link = y~x1+x2+x3)
#' 
#' @export
cancelLP.lvm.reduced  <- function(x, link = NULL, lp = NULL, expar = TRUE, simplify = TRUE, ...){
  x <- lava.reduce:::cancelLP.lvm.reduced(x, link = NULL, lp = NULL, expar = TRUE, restaure = FALSE, simplify = TRUE) 
    
  ## remove penalization
  penalized.link <- penalty(x, type = "link", nuclear = FALSE)
  if(any(penalized.link %in% coef(x) == FALSE)){
    x <- cancelPenalty(x, link = penalized.link[penalized.link %in% coef(x) == FALSE], simplify = simplify)
  }
  penalized.link <- penalty(x, type = "link", nuclear = TRUE)
  if(any(penalized.link %in% coef(x) == FALSE)){
    x <- cancelPenalty(x, link = penalized.link[penalized.link %in% coef(x) == FALSE], simplify = simplify)
  }
  
  return(x)
}

