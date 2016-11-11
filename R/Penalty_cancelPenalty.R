#' @title Remove penalty from a penalized latent variable model
#' @name cancelPenalty
#' @description Remove one or several penalties from a penalized latent variable model
#' 
#' @param x \code{plvm}-object
#' @param link the penalty that should be removed
#' @param value the penalty that should be removed
#' @param simplify if the object contain no more penalty should it be converted to a non penalized lvm
#' 
#' @examples 
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m)
#' 
#' cancelPenalty(pm, link = "y~x5")
#' cancelPenalty(pm) <- "y~x1"
#' cancelPenalty(pm) <- y~x1
#' @export
`cancelPenalty` <-
  function(x,...) UseMethod("cancelPenalty")

#' @export
`cancelPenalty<-` <- function (x, ..., value) {
  UseMethod("cancelPenalty<-", x)
}

`cancelPenalty.plvm` <- function(x, simplify = TRUE, link){
  cancelPenalty(x, simplify) <- link
  return(x)
}

`cancelPenalty<-.plvm` <- function(x, simplify = TRUE, value){
  
  penalty <- penalty(x, type = NULL)
  cancelPenalty(penalty) <- value
  x$penalty <- penalty
  
  if(simplify && length(penalty(x, type = "link")) == 0){
    class(x) <- setdiff(class(x), "plvm")
  }
  
  return(x)
  
}

`cancelPenalty<-.penaltyL12` <- function(x, value){
  
  link <- x$link
  index.rm <- which(link %in% value)
  
  x$group <- x$group[-index.rm]
  
  x$link <- setdiff(link, value)
  x$V[value,] <- 0
  x$V[,value] <- 0
  
  return(x)
}

# `cancelPenalty<-.penaltyNuclear` TODO!!!!

