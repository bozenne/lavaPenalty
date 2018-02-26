### Penalty_clean.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar 10 2017 (13:23) 
## Version: 
## last-updated: mar 14 2017 (17:39) 
##           By: Brice Ozenne
##     Update #: 24
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ doc
#' @title Simplify a lvm object
#' @name clean.penalty
#' @description Remove variables with no link and simplify the class of the lvm object
#' 
#' @param x \code{lvm.penalized}-object
#' @param simplify.penalty  should the class \code{lvm.penalized} be droped if there is no penalty in the object?
#' @param simplify  should the class of the object be simplified ? Overwrite the simplify.x arguments.
#' @param ... additional arguments to lower level functions
#'
#' @details
#' simplify means remove the class \code{"lava.penalized"} if no linear predictor is in the object.
#' 
#' @examples
#'
#' m <- lvm(Y ~ X1 + X2)
#' pm <- penalize(m)
#' pm
#'
#' cancelPenalty(pm) <- Y~X1
#' cancelPenalty(pm) <- Y~X2
#' pm
#'
#' clean(pm)
#' 
# }}}

#' @rdname clean.penalty
#' @export
lavaPenalty.clean.hook <- function(x, simplify.penalty = TRUE, simplify, ...){

    if("plvm" %in% class(x)){
    
        if(!missing(simplify)){
            simplify.penalty <- simplify
        }

        test.penaltyL12 <- NROW(penalty(x, type = "Vlasso")$Vlasso)+ NROW(penalty(x, type = "Vridge")$Vridge)>0
        test.penaltyNuclear <- length(penalty(x, type = "link", nuclear = TRUE))>0
    
        if(simplify.penalty && test.penaltyL12 == FALSE && test.penaltyNuclear == FALSE){
            x$penalty <- NULL
            x$penaltyNuclear <- NULL
            class(x) <- setdiff(class(x), "plvm")
        }
    }

    return(list(x=x))
}




#----------------------------------------------------------------------
### Penalty_clean.R ends here
