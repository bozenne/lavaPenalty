#' @export
`coef0` <-
  function(x,...) UseMethod("coef0")

#' @title Extract Model Coefficients
#'
#' @description Extract the 0 or non 0 coefficients
#'
#' @param x a penalized lvm model
#' @param tol the threshold below which a coefficient is considered to be null
#' @param operator the operator used to extract the coefficients
#' @param penalized should only the penalized coefficient be returned
#' @param value should the value of the coefficients be returned? Otherwise their position.
#' 
#' @examples 
#' set.seed(10)
#' n <- 300
#' formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:12), collapse = "+")))
#' mSim <- lvm(formula.lvm)
#' df.data <- sim(mSim,n)
#' 
#' pm <- penalize(mSim, c("Y~X1","Y~X4", "Y~X10"))
#' pm.fit <- estimate(pm, lambda1 = 1e5, data = df.data, control = list(constrain = TRUE))
#' coef0(pm.fit)
#' coef0(pm.fit, operator = ">")
#' coef0(pm.fit, operator = ">", penalized = TRUE)
#' 
#' @export
coef0.plvmfit <- function(x, tol = 0, operator = "<=", penalized = FALSE, value = TRUE){
  
  names.coef <- names(coef(x)) 
  
  if(!is.null(tol)){
    coefTempo <- names.coef[do.call(operator, args = list(abs(coef(x)), tol))]
  }else{
    coefTempo <- names.coef
  }
  
  if(penalized){
    coefTempo <- intersect(coefTempo, x$penalty$name.coef)
  }
  
  if(value){
    return(coef(x)[coefTempo])
  }else{
    return(coefTempo)
  }
}
