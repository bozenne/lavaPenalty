`coefVar` <-
  function(x,...) UseMethod("coefVar")

`coefCov` <-
  function(x,...) UseMethod("coefCov")

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
  names.var <- paste(x$index$vars, x$index$vars, sep = ",")
  return(grep(paste(names.var, collapse = "|"), coef(x), value = value))
}

#' @title Extract the name or the position of the covariance coefficients
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
#' coefCov(m)
#' coefCov(m, value = TRUE)
#' 
#' covariance(m) <- y1 ~ y2
#' coefCov(m)
#' coefCov(m, value = TRUE)
#' 
#' @export
coefCov.lvm <- function(x, value = FALSE){
  names.cov <- setdiff(coef(x)[x$index$parBelongsTo$cov],
                       paste(x$index$vars, x$index$vars, sep = ",")
  )
  
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
  
  expr <- paste("^",endogenous(x),"\\~", collapse = "|", sep = "")
  index.loadings <- grep(expr, names(coef(x)), value = TRUE)
  loadings <- summary(x)$coef[index.loadings,]
  
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