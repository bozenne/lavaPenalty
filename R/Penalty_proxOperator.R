
#' @title Proximal operators
#' @name proximalOperators
#' @aliases prox proxL1 proxL2 proxE2 proxNuclear
#' 
#' @description Proximal operators for the lasso, ridge, group lasso and nuclear norm penalties.
#' 
#' @param x value of the coefficients
#' @param step value of the step size. Should be between 0 and the inverse of the Lipschitw constant of grad f
#' @param lambda value of the penalization parameter
#' @param nrow the number of rows of the matrix of coefficients
#' @param ncol the number of columns of the matrix of coefficients
#'
#' @return the updated vector of coefficients

#' @rdname proximalOperators

# {{{ proxL1: lasso
proxL1 <- function(x, step, lambda){
    max(0, (abs(x) - lambda * step)) * sign(x)
    # valid result even when lambda is null
}
# }}}

# {{{ proxL2: ridge
#' @rdname proximalOperators
proxL2 <- function(x, step, lambda){
    return( 1 / (1 + lambda * step) * x )
    # valid result even when lambda is null
}
# }}}

# {{{ proxE2: group lasso
#' @rdname proximalOperators
proxE2 <- function(x, step, lambda){
    return( max(0, 1 - sqrt(length(x)) * lambda * step/norm(x, type = "2"))*x )
    # valid result even when lambda is null
    # compared to the documation there is a factor sqrt(length(x))
    # to have a fair penalty for groups of variables with different size
}
# }}}

# {{{  proxNuclear: nuclear norm
#' @rdname proximalOperators
proxNuclear <-  function(x, step, lambda, nrow, ncol){  
  eigen.Mx <- svd(matrix(x, nrow = nrow, ncol = ncol))
  n.eigen <- min(nrow, ncol)
  b <- mapply(proxL1, x = eigen.Mx$d, step = step, lambda = rep(lambda, n.eigen), test.penalty = rep(1, n.eigen))
  return(  as.vector(eigen.Mx$u %*% diag(b) %*% t(eigen.Mx$v)) )
}
# }}}


