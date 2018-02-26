### Penalty_kill.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr  5 2017 (14:42) 
## Version: 
## last-updated: 
##           By: Brice Ozenne
##     Update #: 13
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ doc
#' @title Remove penalty after removing a variable
#' @description Remove penalty after removing a variable
#' @name kill.penalty
#' 
#' @param x \code{plvm}-object
#' @param var the names of the variables that should be removed
#' @param value the names of the link that should be removed
#' @param ... not used
#'
#' @examples
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m)
#' kill(pm) <- "x7"
#' 
#' res <- cancel(pm, y~x1+x2)
#' res
#' 
#' res <- cancel(pm, y~x1+x2+x3)
#' res

# }}}

# {{{ lavaPenalty.remove.hook
#' @rdname kill.penalty
#' @export
lavaPenalty.remove.hook  <- function(x, var, ...){

    if("plvm" %in% class(x)){
        dt.penalty <- penalty(x)
        index.rm <- union(which(dt.penalty$endogenous %in% var),
                          which(dt.penalty$exogenous %in% var))
        if(length(index.rm)>0){
        penalty.rm <- unique(dt.penalty$link[index.rm])
        cancelPenalty(x) <- penalty.rm
        }
    }
    return(x)
}
# }}}

# {{{ lavaPenalty.cancel.hook
lavaPenalty.cancel.hook <- function(x, ..., value){

    if("plvm" %in% class(x)){

        ## normalize value        
        if(is.character(value) && all(value %in% vars(x,lp=TRUE,xlp=TRUE))){
            ## arguments coming from cancel.lvm
            allCoef <- initVar_link(value[1],value[2], format = "txt.formula")
        }else{
            ## argument coming from lavaReduce.remove.hook
            allCoef <- initVar_links(value, format = "txt.formula")
        }
    
        ## remove penalties associated to the link
        dt.penalty <- penalty(x)
        if(any(allCoef %in% dt.penalty$link)){
            cancelPenalty(x) <- allCoef[allCoef %in% dt.penalty$link]
        }
   
    }
    return(x)
}


# }}}


#----------------------------------------------------------------------
### Penalty_kill.R ends here
