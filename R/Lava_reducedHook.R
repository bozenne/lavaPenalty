#' @title Hook to estimate a reduce lvm model
#' @description Add LP to dataset, update the estimator for handling LP, and find initialisation.
#' 
#' @example 
#' R/examples/EX_BuyseTest.R
lava.reduced.estimate.hook <- function(x,data,weight,weight2,estimator,...) {
  
  dots <- list(...)
  if("lvm.reduced" %in% class(x) && length(lp(x))>0){
    
    ## add observations for Linear predictors (to update data with the LP column)
    name.LP <- lp(x)
    data <- cbind(data,
                  data.frame(matrix(0,nrow = NROW(data), ncol = length(name.LP), dimnames = list(NULL,name.LP)))
    )
    
    ## update estimator
    if(estimator == "gaussian"){
      estimator <- "gaussian1"
    }
    
    validEstimator<- paste0("gaussian",c("",1,2))
    if(estimator %in% validEstimator){
      estimator <- paste0(estimator,"LP")
    } else {
      stop("reduced estimator for ",estimator," not implemented \n",
           "available estimators: ",paste(validEstimator, collapse = " "),"\n")
    }
    
    ## intialisation of the parameters
    test.plvm <- ("plvm" %in% class(x)) && ("penalty" %in% names(dots$optim))
    if(is.null(dots$optim$start) && test.plvm == FALSE){
      
      ## non LP
      x0 <- cancelLP(x, simplify = TRUE) # remove LP and external parameters
      startLVM <- coef(estimate(x0, data)) # starting value for the submodel without lp (i.e. values set to 0)
      
      ## LP
      # ls.coef <- lapply(names(x$lp), function(j){ # linear regression to initialize the value of the lp
      #   coef <- lm.fit(y = data[[j]], x = as.matrix(data[x$lp[[j]]$x]))$coefficients
      #   names(coef) <- paste(j,names(coef),sep="~")
      #   return(coef)
      # })
      # startLP <- unlist(ls.coef)[na.omit(match(coef(x),names(unlist(ls.coef))))]
      
      ## update
      dots$optim$start <- setNames(rep(0, length = length(coef(x))), coef(x))
      # dots$optim$start[names(startLP)] <- startLP
      dots$optim$start[names(startLVM)] <- startLVM
      dots$optim$start <- dots$optim$start[which(!is.na(dots$optim$start))]
    }
    
  }
  
  
  return(c(list(x=x,data=data,weight=weight,weight2=weight2,estimator=estimator),dots)) 
}

#' @description add the LPs in the first and second empirical moment so that LP is not converted to a latent variable
procdata.lvm.reduced <- function(x, data, ...){
  
  name.LP <- lp(x)
  data <- cbind(data,
                data.frame(matrix(0,nrow = NROW(data), ncol = length(name.LP), dimnames = list(NULL,name.LP)))
  )

  class(x) <- setdiff(class(x),"lvm.reduced")
  return(procdata(x, data=data, ...))
}


lava.reduced.post.hook <- function(x){
  
  return(x)
  
}