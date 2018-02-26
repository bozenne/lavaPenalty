# {{{ doc
#' @title Remove penalty from a penalized latent variable model
#' @name cancelPenalty
#' @description Remove one or several penalties from a penalized latent variable model
#' 
#' @param x \code{plvm}-object
#' @param value the penalty that should be removed
#' @param rm.lasso should lasso penalties be removed?
#' @param rm.ridge should ridge penalties be removed?
#' @param rm.groupLasso should group lasso penalties be removed?
#' @param extraParameter external parameters
#' @param ... additional arguments to be passed to lower level functions
#'
#' @details
#' Since \code{lavaReduce::initVar_links} does not work for the external parameters, they have to be handled separately.
#' 
#' @examples
#' ## lasso
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m)
#' 
#' cancelPenalty(pm, value = "y~x5")
#' cancelPenalty(pm) <- y~x1
#' cancelPenalty(pm) <- c("y~x2","y~x3")
#' cancelPenalty(pm) <- penalty(pm, no.ridge = TRUE)$link
#' pm
#'
#' ## group lasso
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' categorical(m, K = 3, labels = 1:3) <- ~x1
#' pm <- penalize(m)
#' pm
#' cancelPenalty(pm) <- "y~x2"
#' pm
#' cancelPenalty(pm) <- "y~x12"
#' pm
#' cancelPenalty(pm) <- penalty(pm, no.lasso = TRUE, no.ridge = TRUE)$link
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
`cancelPenalty.plvm` <- function(x, ..., value){
  cancelPenalty(x, ...) <- value
  return(x)
}
# }}}
# {{{ cancelPenalty<-.plvm
#' @rdname cancelPenalty
`cancelPenalty<-.plvm` <- function(x, ..., value){

    penalty <- penalty(x, type = "object")
    cancelPenalty(penalty, extraParameter = coefExtra(x, value = TRUE), ...) <- value
    x$penalty <- penalty
      
  return(x)
  
}
# }}}
# {{{ cancelPenalty<-.penaltyL12
#' @rdname cancelPenalty
`cancelPenalty<-.penaltyL12` <- function(x, extraParameter,
                                         rm.lasso = TRUE, rm.ridge = TRUE, rm.groupLasso = TRUE,
                                         value){
    ## deal with external parameters
    if(is.character(value)){
        value.EP <- intersect(value, extraParameter)
        value <- setdiff(value, extraParameter)
    }else{
        value.EP <- NULL
    }
    
    ## normalize argument value
    # initVar_links cannot deal with external parameters like x1:0|1 since it is not a formula
    value.P <- lavaReduce::initVar_links(value,
                                          format = "txt.formula")    
    value <- c(value.P, value.EP)

    ## identify all penalties
    table.penalty <- penalty(x)
    if(any(value %in% table.penalty$link == FALSE)){
        stop("Cannot remove an non existing link in object \n",
             "link: ",paste(value[value %in% table.penalty$link == FALSE], collapse = " "),"\n")
    }

    if(rm.lasso && table.penalty[link %in% value & penalty == "lasso",.N]>0){
        link.lasso <- table.penalty[link %in% value & penalty == "lasso",link]

        Vtempo <- penalty(x, type = "Vlasso")$Vlasso
        Vtempo[rownames(Vtempo) %in% link.lasso,] <- 0
        indexN0 <- which(Matrix::colSums(abs(Vtempo))!=0)        
        
        lambda1 <- penalty(x, type = "lambda1")$lambda1

        penalty(x, type = "Vlasso", add = FALSE) <- Vtempo[,indexN0, drop = FALSE]
        if(length(lambda1)>0){
            penalty(x, type = "lambda1", add = FALSE) <- lambda1[indexN0]
        }
    }
    if(rm.ridge && table.penalty[link %in% value & penalty == "ridge",.N]>0){
        link.ridge <- table.penalty[link %in% value & penalty == "ridge",link]
        
        Vtempo <- penalty(x, type = "Vridge")$Vridge
        Vtempo[rownames(Vtempo) %in% link.ridge,] <- 0
        indexN0 <- which(Matrix::colSums(abs(Vtempo))!=0)

        lambda2 <- penalty(x, type = "lambda2")$lambda2
        
        penalty(x, type = "Vridge", add = FALSE) <- Vtempo[,indexN0, drop = FALSE]
         if(length(lambda2)>0){
            penalty(x, type = "lambda2", add = FALSE) <- lambda2[indexN0]
        }
    }
    if(rm.groupLasso && table.penalty[link %in% value & penalty == "group lasso",.N]>0){
        link.group <- table.penalty[link %in% value & penalty == "group lasso",link]
        
        Vtempo <- penalty(x, type = "Vgroup")$Vgroup
        Vtempo[rownames(Vtempo) %in% link.group,] <- 0
        indexN0 <- which(Matrix::colSums(abs(Vtempo))!=0)
        
        lambdaG <- penalty(x, type = "lambdaG")$lambdaG

        penalty(x, type = "Vgroup", add = FALSE) <- Vtempo[,indexN0, drop = FALSE]
        if(length(lambdaG)>0){
            penalty(x, type = "lambdaG", add = FALSE) <- lambdaG[indexN0]
        }
    }
    
  return(x)
}

# }}}

# {{{ cancelPenalty<-.penaltyNuclear [TODO]
# }}}


