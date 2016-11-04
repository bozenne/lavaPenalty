'.onLoad' <- function(lib, pkg="lava.penalty") {
  
  lava::addhook("lava.penalty.estimate.hook", hook = "estimate.hooks")
  lava::addhook("lava.penalty.post.hook", hook = "post.hooks")
  
  lava::lava.options(proxGrad = list(method = "ISTA", iter.max = 1000, step = 1, BT.n = 100, BT.eta = 0.8, abs.tol = 1e-9, rel.tol = 1e-10, force.descent = FALSE, export.iter = FALSE),
                     EPSODE = list(resolution_lambda1 = c(1e-1,1e-3), nstep_max = min(length(beta)*50,1e4), ode.method = "euler", tol.0 = 1e-8, reversible = FALSE, exportAllPath = FALSE),
                     calcLambda = list(fit = "BIC", warmUp = FALSE), 
                     constrain = TRUE)
}

'.onAttach' <- function(lib, pkg="lava.penalty") {
  desc <- utils::packageDescription(pkg)
  packageStartupMessage(desc$Package, " version ",desc$Version)
}