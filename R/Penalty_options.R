### Penalty_options.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar  7 2017 (15:51) 
## Version: 
## last-updated: apr  4 2017 (12:58) 
##           By: Brice Ozenne
##     Update #: 29
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Set global options for lavaPenalty
#' 
#' @description Extract and set global parameters of lavaPenalty.
#' 
#' @param ... Arguments
#' @param type which set of parameters to update: \code{"lava"}, \code{"proxGrad"}, or \code{"EPSODE"}.
#'
#' @examples
#'
#' ## lava
#' lavaPenalty.options()
#' lavaPenalty.options(trace = TRUE)
#'
#' ## lavaPenalty
#' lavaPenalty.options(type = "proxGrad")
#' lavaPenalty.options(trace = 2, type = "proxGrad")
#' @export
lavaPenalty.options <- function(..., type = "lava"){
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
