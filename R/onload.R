#' @title Default Options for lavaPenalty
#' @description Default options for lavaPenalty.
#' @name onload
#'
#' @details Add the following elements to \code{lava.options}:
#' \itemize{
#' item iter.max maximum number of iterations
#' item trace should the convergence diagnostics be displayed at each step
#' item abs.tol convergence is the difference in likelihood between two consecutive steps is below this threshold
#' item rel.tol convergence is the relative difference  in likelihood between two consecutive steps is below this threshold
#' item step maximun step for the proximal gradient algorithm. 
#' If NULL the step is estimated using the inverse of the maximal eigenvalue of the hessian (in absolute value) and re-estimated at each step
#' Otherwise backtracking is used.
#' item BT.n number of backtracking steps
#' item BT.eta multiplicative factor for the step 
#' item force.descent 
#' item export.iter should all iterations be exported (details.cv)
#' item method type of iteration
#' ISTA correspond to the ISTA step as described in Bech 2009
#' FISTA_Beck correspond to the FISTA step as described in Bech 2009
#' FISTA_Vand correspond to ??
#' mFISTA_Vand correspond to ??
#' item resolution_lambda1 the first value is the maximum relative difference in parameter between two steps. 
#' If this lead to a too small step, the second value is used as the minimum change in penalization parameter between two steps.
#' item increasing direction of the path
#' item stopLambda if not null, stop the path when the penalty parameter reach this value
#' item stopParam if not null, stop the path when the number of 0 (if increasing = TRUE) or non 0 (if increasing = FALSE) parameters has reached this value.
#' item nstep_max the maximum number of iterations
#' item ode.method the type of method to use to solve the ode (see the documentation of deSolve:::ode)
#' item constrain the constrain on the variance parameters
#' item reversible should the algorithm allow a 0 parameter to become non-0 (when increasing = TRUE) or non-0 parameter to become 0 (when increasing = FALSE)
#' item tol.0 tolerance for classify a parameter from beta_lambdaMax in the set of 0 parameters
#' item exportAllPath export all the regularization path (and not only the breakpoints)
#' item trace should the function be traced
#' }
#'

#' @rdname onload
.onLoad <- function(lib, pkg="lavaPenalty") {
  
    lava::addhook("lavaPenalty.estimate.hook", hook = "estimate.hooks")
    lava::addhook("lavaPenalty.post.hook", hook = "post.hooks")
    lava::addhook("lavaPenalty.multigroup.hook", hook = "multigroup.hooks")
    lava::addhook("lavaPenalty.remove.hook", hook = "remove.hooks")
    lava::addhook("lavaPenalty.cancel.hook", hook = "cancel.hooks")
    
    lava::lava.options(proxGrad = list(method = "ISTA", iter.max = 1000, step = 1, BT.n = 100, BT.eta = 0.8, abs.tol = 1e-9, rel.tol = 1e-10, force.descent = FALSE,
                                       export.iter = FALSE, trace = 1),
                       EPSODE = list(resolution_lambda1 = c(1e-1,1e-3), nstep_max = min(length(beta)*50,1e4), ode.method = "euler", tol.0 = 1e-8,
                                     lars = FALSE, reversible = FALSE, increasing = FALSE, stopLambda = NULL, stopParam = NULL,
                                     exportAllPath = FALSE, trace = 1),
                       proxGradPath = list(warmUp = FALSE),
                       calcLambda = list(fit = "BIC"),
                       Nuclear = list(symbols = c("Image[","]")),
                       constrain = TRUE)

    lavaReduce::addhook_lavaReduce("lavaPenalty.clean.hook", hook = "clean.hooks")
    lavaReduce::addhook_lavaReduce("lavaPenalty.reduce.hook", hook = "reduce.hooks")
}

#' @rdname onload
.onAttach <- function(lib, pkg="lavaPenalty") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}

