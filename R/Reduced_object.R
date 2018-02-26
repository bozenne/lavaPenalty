# {{{ doc
#' @title Aggregate exogeneous variables associated with penalties
#' @description Aggregate exogeneous variables associated with penalties into a linear predictor
#' @name reduce.penalty
#' 
#' @param x \code{lvm}-object
#' @param link the links that should be included in the linear predictor
#' @param ... other arguments to be passed to \code{lavaReduce::reduce.lvm}
#' 
#' @examples  
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m)
#' res <- reduce(pm)
#' reduce(pm, link = y~x1+x2+x3)
#' 
# }}}

# {{{ lavaPenalty.reduce.hook
#' @rdname reduce
#' @export
lavaPenalty.reduce.hook <- function(x, link = NULL, ...){

    if("plvm" %in% class(x)){
        if(is.null(link)){
            link <- c(unique(penalty(x, nuclear = FALSE)$link),
                      penalty(x, type = "link", nuclear = TRUE)
                      )
        }
        save <- x[c("penalty","penaltyNuclear")]
    }
    
    return(list(link = link,
                save = save))
}
# }}}

