#' @param iter.max maximum number of iterations
#' @param trace should the convergence diagnostics be displayed at each step
#' @param abs.tol convergence is the difference in likelihood between two consecutive steps is below this threshold
#' @param rel.tol convergence is the relative difference  in likelihood between two consecutive steps is below this threshold
#' @param step maximun step for the proximal gradient algorithm. 
#' If NULL the step is estimated using the inverse of the maximal eigenvalue of the hessian (in absolute value) and re-estimated at each step
#' Otherwise backtracking is used.
#' @param BT.n number of backtracking steps
#' @param BT.eta multiplicative factor for the step 
#' @param force.descent 
#' @param export.iter should all iterations be exported (details.cv)
#' @param method type of iteration
#' ISTA correspond to the ISTA step as described in Bech 2009
#' FISTA_Beck correspond to the FISTA step as described in Bech 2009
#' FISTA_Vand correspond to ??
#' mFISTA_Vand correspond to ??
#' @param resolution_lambda1 the first value is the maximum relative difference in parameter between two steps. 
#' If this lead to a too small step, the second value is used as the minimum change in penalization parameter between two steps.
#' @param increasing direction of the path
#' @param stopLambda if not null, stop the path when the penalty parameter reach this value
#' @param stopParam if not null, stop the path when the number of 0 (if increasing = TRUE) or non 0 (if increasing = FALSE) parameters has reached this value.
#' @param nstep_max the maximum number of iterations
#' @param ode.method the type of method to use to solve the ode (see the documentation of deSolve:::ode)
#' @param constrain the constrain on the variance parameters
#' @param reversible should the algorithm allow a 0 parameter to become non-0 (when increasing = TRUE) or non-0 parameter to become 0 (when increasing = FALSE)
#' @param tol.0 tolerance for classify a parameter from beta_lambdaMax in the set of 0 parameters
#' @param exportAllPath export all the regularization path (and not only the breakpoints)
#' @param trace should the function be traced
#' 

.onLoad <- function(lib, pkg="lava.penalty") {
  
  lava::addhook("lava.penalty.estimate.hook", hook = "estimate.hooks")
  lava::addhook("lava.penalty.post.hook", hook = "post.hooks")
  
    lava::lava.options(proxGrad = list(method = "ISTA", iter.max = 1000, step = 1, BT.n = 100, BT.eta = 0.8, abs.tol = 1e-9, rel.tol = 1e-10, force.descent = FALSE,
                                       export.iter = FALSE, trace = 1),
                       EPSODE = list(resolution_lambda1 = c(1e-1,1e-3), nstep_max = min(length(beta)*50,1e4), ode.method = "euler", tol.0 = 1e-8,
                                     reversible = FALSE, increasing = FALSE, stopLambda = NULL, stopParam = NULL,
                                     exportAllPath = FALSE, trace = 1),
                     calcLambda = list(fit = "BIC", warmUp = FALSE),
                     Nuclear = list(symbols = c("Image[","]")),
                     constrain = TRUE)  
}

.onAttach <- function(lib, pkg="lava.penalty") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}

lava_matrices.lvm <- get("matrices.lvm", envir = asNamespace("lava"), inherits = FALSE)
lava_categorical2dummy <- get("categorical2dummy", envir = asNamespace("lava"), inherits = FALSE)
lava_estimate.lvm <- get("estimate.lvm", envir = asNamespace("lava"), inherits = FALSE)
