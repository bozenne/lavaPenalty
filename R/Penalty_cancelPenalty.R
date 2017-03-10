# {{{ cancelPenalty

# {{{ doc
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
#'
#' ## lasso
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m)
#' 
#' cancelPenalty(pm, link = "y~x5")
#' cancelPenalty(pm) <- y~x1
#' cancelPenalty(pm) <- c("y~x2","y~x3")
#' cancelPenalty(pm) <- penalty(pm)
#' pm
#'
#' ## group lasso
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' categorical(m, K = 3) <- ~x1
#' pm <- penalize(m)
#' cancelPenalty(pm) <- "x1:0|1"
#' pm
#' cancelPenalty(pm) <- "x1:1|2"
#' pm
#' cancelPenalty(pm) <- penalty(pm, type = "link")
#' pm
#' 
#' 
#' @export
`cancelPenalty` <-
    function(x,...) UseMethod("cancelPenalty")
# }}}
# {{{ cancelPenalty<-
#' @rdname cancelPenalty
#' @export
`cancelPenalty<-` <- function (x, ..., value) {
  UseMethod("cancelPenalty<-", x)
}
# }}}
# {{{ cancelPenalty.plvm
#' @rdname cancelPenalty
`cancelPenalty.plvm` <- function(x, simplify = TRUE, link){
  cancelPenalty(x, simplify) <- link
  return(x)
}
# }}}
# {{{ cancelPenalty<-.plvm
#' @rdname cancelPenalty
`cancelPenalty<-.plvm` <- function(x, clean = TRUE, ..., value){
  
    penalty <- penalty(x, type = NULL)
    cancelPenalty(penalty, extraParameter = coefExtra(pm, value = TRUE)) <- value
    x$penalty <- penalty
    if(clean){
        x <- clean(x, ...)
    }
  
  return(x)
  
}
# }}}
# {{{ cancelPenalty<-.penaltyL12
#' @rdname cancelPenalty
`cancelPenalty<-.penaltyL12` <- function(x, extraParameter, value){

    ## normalize argument value
    # initVar_links cannot deal with external parameters like x1:0|1 since it is not a formula
    value.P <- lava.reduce::initVar_links(setdiff(value, extraParameter),
                                          format = "txt.formula")
    value.EP <- intersect(value, extraParameter)
    value <- c(value.P, value.EP)

    ##
    link.elasticNet <- penalty(x, type = "link", no.group = TRUE)
    if(!is.null(link.elasticNet)){        
        value.elasticNet <- value[value %in% link.elasticNet]
    }else{
        value.elasticNet <- NULL
    }
    
    link.groupLasso <- penalty(x, type = "link", no.elasticNet = TRUE)
    if(!is.null(link.groupLasso)){
        value.groupLasso <- value[value %in% link.groupLasso]
    }else{
        value.groupLasso <- NULL
    }

    
    if(length(value) != length(value.elasticNet)+length(value.groupLasso)){
        wrongValue <- value[value %in% c(value.elasticNet,value.groupLasso) == FALSE]
        stop("Cannot remove an non existing link in object \n",
             "link: ",paste(wrongValue, collapse = " "),"\n")
    }

    if(length(value.elasticNet)>0){        
        Vtempo <- penalty(x, type = "VelasticNet")
        Vtempo[rownames(Vtempo) %in% value.elasticNet,] <- 0
        penalty(x, type = "VelasticNet") <- Vtempo[,Matrix::colSums(abs(Vtempo))!=0, drop = FALSE]
    }
    
    if(length(value.groupLasso)>0){
        Vtempo <- penalty(x, type = "Vgroup")
        Vtempo[rownames(Vtempo) %in% value.groupLasso,] <- 0
        penalty(x, type = "Vgroup") <- Vtempo[,Matrix::colSums(abs(Vtempo))!=0, drop = FALSE]
    }
    
  return(x)
}
# }}}

# {{{ cancelPenalty<-.penaltyNuclear [TODO]
# }}}

# }}}

