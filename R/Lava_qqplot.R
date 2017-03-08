# {{{ residplot

#' @title Plot the residuals of a lvm object
#' @name residplot
#'
#' @param object a lvm model
#' @param data a data frame
#' @param factor.vcov inflation factor for the variance when sampling the initialization points
#' @param n.init number of initialization points to be used
#' @param ncpus the number of CPU to be used
#' @param keep.cov should the covariance between parameter be kept to simulate the initialization points
#' @param verbose should a progression bar be displayed?
#' @param ... additional arguments to be passed to estimate
#' 
#' @details 
#' Simulation is based on a multivariate truncated normal law (even though it is not satifying for the variance components)
#' 
#' @return a data frame/cvlvm object containing the convergence status (by default 0 indicates successful convergence, see ?optim), the value of the log-likelihood and the estimated parameters (in columns) for each initialization (in rows)
#' 
#' @examples 
#' m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
#' covariance(m) <- v1~v2+v3+v4
#' latent(m) <- ~ x
#' dd <- sim(m,100) ## Simulate 100 observations from model
#' e <- estimate(m, dd) ## Estimate parameters
#'
#' residplot(e)
#' 
#' @export
residplot <- function (x, ...) {
  UseMethod("residplot", x)
}

#' @rdname cvCheck
#' @export
residplot.lvmfit <- function(object, variables = NULL, mfrow = NULL,
                             type = "qqtest",  centralPercents = 0.95,...){

    M.res <- predict(object, residual = TRUE)
    name.vars <- colnames(M.res)
    
    if(!is.null(variables)){
        if(any(variables %in% name.vars == FALSE)){
            stop("unknown variable(s): ",paste(variables[variables %in% name.vars == FALSE], collapse = " "),"\n",
                 "endogenous variables: ",paste(endogenous(object), collapse = " "),"\n",
                 "latent variables: ",paste(latent(object), collapse = " "),"\n")
        }
        M.res <- M.res[,variables,drop=FALSE]
        name.vars <- variables
    }
    if(type %in% c("qqtest","qqnorm") == FALSE){
        stop("wrong specification of type \n",
             "must be \"qqtest\" or \"qqnorm\" \n")
    }

    n.var <- NCOL(M.res)       
    if(is.null(mfrow)){
        mfrow <- c(round(sqrt(n.var)), ceiling(n.var/round(sqrt(n.var))))
    }
    op <- par(mfrow = mfrow)
    sapply(1:n.var, function(row){
        resid <- na.omit(M.res[,row])
        main <- name.vars[row]
        if(all(resid < 1e-5)){
            plot(0,0, col = "white", axes = FALSE, xlab = "", ylab = "", main = main)
            text(0,0,"all residuals < 1e-5")
        }else if(type == "qqtest"){
            qqtest::qqtest(resid, main = name.vars[row],
                           centralPercents = centralPercents,
                           ...)
        }else if(type == "qqnorm"){
            qnorm(M.res[,row], main = name.vars[row])
        }
    })
    par(op)

    return(invisible(M.res))
}



# }}}
