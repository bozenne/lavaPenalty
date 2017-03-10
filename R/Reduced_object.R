# {{{ reduce.plvm
#' @title Reduce latent variable model
#' @name reduce
#' @description Aggregate exogeneous variables into a linear predictor
#' 
#' @param x \code{lvm}-object
#' @param endo the endogeneous variables for which the related exogeneous variables should be aggregated
#' @param rm.exo should the exogeneous variables be remove from the object
#' @param ... other arguments to be passed to \code{lava.reduce::reduce.lvm}
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
reduce.plvm <- function(object, link = NULL, ...){
  
    if(is.null(link)){
        link <- c(penalty(object, type = "link", nuclear = FALSE),
                  penalty(object, type = "link", nuclear = TRUE)
                  )
    }
    object <- lava.reduce:::reduce.lvm(object, link = link, ...)
  
    return(object)
}
# }}}

