#' @export
`coefVar` <-
  function(x,...) UseMethod("coefVar")

#' @export
`coefCov` <-
  function(x,...) UseMethod("coefCov")

#' @export
`loadings` <- function(object, ...) UseMethod("loadings")

#' @title Extract the name or the position of the variance coefficients
#' 
#' @param x a lvm model
#' @param value should the name of the coefficient be returned? Else return the coefficients
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x
#' latent(m) <- ~u
#' 
#' coefVar(m)
#' coefVar(m, value = TRUE)
#' 
#' @export
coefVar.lvm <- function(x, value = FALSE){ 
  names.var <- paste(x$index$vars, x$index$vars, sep = lava.options()$symbol[2])
  return(grep(paste(names.var, collapse = "|"), coef(x), value = value))
}

#' @title Extract the name or the position of the covariance coefficients
#' 
#' @param x a lvm model
#' @param value should the name of the coefficient be returned? Else return the coefficients
#' @param keep.var should the variance parameters be output?
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x
#' latent(m) <- ~u
#' 
#' coefCov(m)
#' coefCov(m, value = TRUE)
#'
#' coefCov(m, keep.var = FALSE)
#' coefCov(m, value = TRUE, keep.var = FALSE)
#' 
#' covariance(m) <- y1 ~ y2
#' coefCov(m)
#' coefCov(m, value = TRUE)
#' 
#' coefCov(m, keep.var = FALSE)
#' coefCov(m, value = TRUE, keep.var = FALSE)
#'
#' 
#' @export
coefCov.lvm <- function(x, value = FALSE, keep.var = TRUE){
  names.cov <- coef(x)[x$index$parBelongsTo$cov]
  if(keep.var == FALSE){
    names.cov <- setdiff(names.cov,
                         coefVar(x, value = TRUE))
  }
  
  if(length(names.cov)==0){
    return(NULL)
  }else{
    return(grep(paste(names.cov, collapse = "|"), coef(x), value = value))
  }
}



#' @title Extract the summary table for the loadings
#' @name loadings
#' @aliases loadings loadings.lvmfit loadings.lvm.missing
#' 
#' @param x a lvm model
#' @param value should the name of the coefficient be returned? Else return the coefficients
#'
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x
#' latent(m) <- ~u
#' df <- sim(m,1e2)
#' 
#' fit <- estimate(m, df)
#' 
#' loadings(fit) 
#' loadings(fit, col = "P-value") 

#' @rdname loadings
#' @export
loadings.lvmfit <- function(x, col = NULL){
  
  # candidates
  index.reg <- x$model$index$parBelongsTo$reg
  
  # check that the first variable is endogeneous and the second a latent variable
  test <- sapply(names(coef(x))[index.reg], function(var){
    vars <- initVar_link(var)
    test <- (vars[1] %in% endogenous(x))*(vars[2] %in% latent(x))
    return(test)
  })
  names.loadings <- names(coef(x)[index.reg[test==1]])
  
  # extraction
  loadings <- summary(x)$coef[names.loadings,]
  
  if(is.null(col)){
    return(loadings)
  }else{
    match.arg(col, choices = c("Estimate","Std. Error", "Z-value","P-value"))
    return(loadings[,col])
  }
}

#' @rdname loadings
#' @export
loadings.lvm.missing <- loadings.lvmfit