#' @title Hook to estimate a reduce lvm model
#' @description Add LP to dataset, update the estimator for handling LP, and find initialisation.
#' 
#' @example 
#' R/examples/EX_reducedModel.R
#' @export
lava.reduced.estimate.hook <- function(x,data,weight,weight2,estimator,...) {
  
  dots <- list(...)
  if("lvm.reduced" %in% class(x) && length(lp(x))>0){
    
    ## set reference endogeneous variable where there is no linear predictor
    if(length(latent(x))>0){
      n.latent <- length(latent(x))
      for(iterLatent in 1:n.latent){
        endoLatent <- names(which(x$index$M[latent(x)[iterLatent],]==1))
        reference <- setdiff(endoLatent, names(x$lp))
        if(length(reference)==0){
          stop("A reduce LVM must have for each latent variable at least one observed outcome not in reduced form \n")
        }else{
          f <- as.formula(paste(reference, latent(x)[iterLatent], sep = "~"))
          regression(x, f) <- 1
          f <- as.formula(paste(reference, 1, sep = "~"))
          intercept(x, f) <- 0
        }
      }
    }
    
    ## add observations for Linear predictors (to update data with the LP column)
    name.LP <- lp(x)
    data <- cbind(data,
                  data.frame(matrix(0,nrow = NROW(data), ncol = length(name.LP), dimnames = list(NULL,name.LP)))
    )
    
    ## 
    test.plvm <- ("plvm" %in% class(x)) && ("penalty" %in% names(dots$optim))
    if(test.plvm){
      
      lambdaN <- penalty(x, nuclear = TRUE, "lambdaN")
      if(!is.null(lambdaN) && lambdaN>0 && estimator == "gaussian"){ # set gaussian2 to simplify the computation of the hessian
        estimator <- "gaussian2"
      }
      
    }else if(is.null(dots$optim$start)){ # intialisation of the parameters
      
      ## non LP
      x0 <- cancelLP(x, simplify = TRUE) # remove LP and external parameters
      startLVM <- estimate(x0, data = data, quick = TRUE) # starting value for the submodel without lp (i.e. values set to 0)
      
      ## update
      dots$optim$start <- setNames(rep(0, length = length(coef(x))), coef(x))
      # dots$optim$start[names(startLP)] <- startLP
      dots$optim$start[names(startLVM)] <- startLVM
      dots$optim$start <- dots$optim$start[which(!is.na(dots$optim$start))]
    }
    
    ## update estimator
    if(estimator == "gaussian"){
      estimator <- "gaussian1"
    }
    
    validEstimator <- paste0("gaussian",c("",1,2))
    if(estimator %in% validEstimator){
      estimator <- paste0(estimator,"LP")
    } else {
      stop("reduced estimator for ",estimator," not implemented \n",
           "available estimators: ",paste(validEstimator, collapse = " "),"\n")
    }
    
    
  }
  
  return(c(list(x=x,data=data,weight=weight,weight2=weight2,estimator=estimator),dots)) 
}

#' @title Preparation of data
#'
#' @description add the LPs in the first and second empirical moment so that LP is not converted to a latent variable
procdata.lvm.reduced <- function(x, data, ...){
  
  name.LP <- lp(x)
  data <- cbind(data,
                data.frame(matrix(0,nrow = NROW(data), ncol = length(name.LP), dimnames = list(NULL,name.LP)))
  )
  
  class(x) <- setdiff(class(x),"lvm.reduced")
  return(procdata(x, data=data, ...))
}


#' @export
lava.reduced.post.hook <- function(x){
  
  return(x)
  
}