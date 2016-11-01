#' @title Add regression association to a latent variable model with linear predictor
#' @description Same as regression for latent variable model with an argument reduce to introduce linear predictors
#' 
#' @param x \code{lvm}-object
#' @param to Character vector of outcome(s) or formula object.
#' @param from Character vector of predictor(s).
#' @param y	Alias for 'to'
#' @param x	Alias for 'from'
#' @param reduce should the from variable be grouped into a linear predictor
#' 
#' @examples 
#' 
#' m <- lvm.reduced()
#' m <- regression(m,y='y1',x='x'%++%1:2)
#' m <- regression(m,y='y1',x='z'%++%1:5, reduce = TRUE)
#' 


regression.lvm.reduced <- function(object = lvm(), to, from, y, x, reduce = FALSE, ...){
  
  if(identical(reduce,TRUE)  || is.character(reduce)){
    
    #### extract to and from
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
    
    #### specific case - see regression.lvm
    if (missing(to)) {
      return(regfix(object))
    }
    if (inherits(to, "formula")) {
      if (!missing(value)) {
        regression(object, to, reduce = reduce, ...) <- value
      }
      else {
        regression(object, reduce = reduce, ...) <- to
      }
      object$parpos <- NULL
      return(object)
    }
    if (is.list(to)) {
      for (t in to) regression(object, reduce = reduce, ...) <- t
      object$parpos <- NULL
      return(object)
    }
    
    #### reduce
    allCoef <- coef(regression(object, to = to, from = from, reduce = FALSE, ...))
    if(!identical(reduce,TRUE) && length(reduce)!=length(to)){
      stop("wrong specification of argument \'reduce\' \n",
           "must be TRUE or a character vector of size ",length(to),"\n")
    }
    
    if("lp" %in% names(object) == FALSE){object$lp <- list()}
    for(iterR in to){ ### I don't know where to get the constrains
      name.LP <- if(reduce==TRUE){paste0("LP",iterR)}else{reduce}
      regression(object, to = iterR, from = name.LP, reduce = FALSE, ...) <- 1
      if(iterR %in% names(object$lp)){ 
        # object$lp[[iterR]]$indexCoef <- match(object$lp[[iterR]]$indexCoef, coef(object))
        object$lp[[iterR]]$x <- c(object$lp[[iterR]]$x, from)
        object$lp[[iterR]]$con <- c(object$lp[[iterR]]$con, rep(NA,length(from)))
      }else{
        # object$lp[[iterR]]$indexCoef <- match(namesCoef, coef(object))
        object$lp[[iterR]]$x <- from
        object$lp[[iterR]]$con <- rep(NA,length(from))
        object$lp[[iterR]]$name <- name.LP
      }
    }
    parameter(object) <- setdiff(allCoef,coef(object))
    return(object)
    
  }else {
    
    return(lava:::regression.lvm(object = object, to = to, from = from, y = y, x = x, ...))
    
  }
  
  
}