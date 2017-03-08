### Penalty_options.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar  7 2017 (15:51) 
## Version: 
## last-updated: mar  7 2017 (16:25) 
##           By: Brice Ozenne
##     Update #: 27
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Set global options for lava.penalty
#' @describe Extract and set global parameters of lava.penalty.
#' 
#' @param ... Arguments
#' @param type which set of parameters to update: \code{"lava"}, \code{"proxGrad"}, or \code{"EPSODE"}.
#'
#' @examples
#'
#' ## lava
#' lava.penalty.options()
#' lava.penalty.options(trace = TRUE)
#'
#' ## lava.penalty
#' lava.penalty.options(type = "proxGrad")
#' lava.penalty.options(trace = 2, type = "proxGrad")
#' @export
lava.penalty.options <- function(..., type = "lava"){
    n.args <- length(list(...))
    additional.types <- c("lava", "proxGrad", "EPSODE", "Nuclear", "calcLambda")
    
    type <- match.arg(type, additional.types)
    if(type == "lava"){
        res <- lava::lava.options(...)
        res[additional.types] <- NULL
    }else{
        if(n.args==0){
            res <- lava::lava.options()[[type]]
        }else{
            res <- lava::lava.options()[[type]]
            dots <- list(...)
            names.dots <- names(dots)
            if(any(names.dots %in% names(res) == FALSE)){
                stop("unknown options for type ",type,"\n",
                     "unknown options: ",paste(names.dots[names.dots %in% names(res) == FALSE]),"\n") 
            }
            res[names.dots] <- dots
            switch(type,
                   "proxGrad" = lava::lava.options(proxGrad = res),
                   "EPSODE" = lava::lava.options(EPSODE = res),
                   "Nuclear" = lava::lava.options(Nuclear = res),
                   "calcLambda" = lava::lava.options(calcLambda = res),
                   )
        }

    }

    if(n.args>0){
        return(invisible(res))
    }else{
        return(res)
    }
}

#----------------------------------------------------------------------
### Penalty_options.R ends here
