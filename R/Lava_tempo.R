checkMoment <- function(lvmRed, elvm, param = coef(e), data = elvm$data$model.frame){
  
  indexRED <- na.omit(match(coef(lvmRed),names(param)))
  
  ## test logLik
  expect_equal(gaussianReduced_logLik.lvm(lvmRed, p = param[indexRED], data = data),
               lava:::gaussian_logLik.lvm(object = elvm, data=data, p = param) 
  )
  
  ## test gradient
  expect_equal(unname(gaussianReduced_gradient.lvm(lvmRed, p = param[indexRED], data = data, indiv = FALSE)),
               lava:::gaussian_gradient.lvm(x = elvm$model, data= data, p=param, S = elvm$S, n = elvm$data$n, mu = elvm$mu)[indexRED]
  )
  
  ## test score
  expect_equal(unname(gaussianReduced_score.lvm(lvmRed, p = param[indexRED], data = data, indiv = FALSE)),
               unname(lava:::gaussian_score.lvm(x = elvm$model, data=data, p=param, S = elvm$S, n = elvm$data$n, mu = elvm$mu)[,indexRED,drop=FALSE])
  )  
  expect_equal(unname(gaussianReduced_score.lvm(lvmRed, p = param[indexRED], data = data, indiv = TRUE)),
               unname(score(elvm$model,data=data,p=param,indiv=TRUE)[,indexRED,drop=FALSE])
  )
  expect_equal(unname(gaussianReduced_score.lvm(lvmRed, p = param[indexRED], data = data, indiv = FALSE)),
               unname(score(elvm$model,data=data, p = param,indiv=FALSE)[,indexRED,drop=FALSE])
  )  
  
  ## test hessian
  Hred_E <- gaussianReduced2_hessian.lvm(x = lvmRed, data=data, p=param[indexRED])
  Hlava <- lava:::gaussian2_hessian.lvm(x = elvm$model, data=data, p=param, n = elvm$data$n, mu = elvm$mu, S = elvm$S)
  expect_equal(attr(Hred_E,"grad"),
               attr(Hlava,"grad")[indexRED]
  ) 
  attr(Hred_E,"grad") <- NULL
  attr(Hlava,"grad") <- NULL
  expect_equal(unname(Hred_E), unname(Hlava[indexRED,indexRED]))
  
  Hred_num <- gaussianReduced1_hessian.lvm(x = lvmRed, data=data, p=param[indexRED])
  Hlava <- lava:::gaussian1_hessian.lvm(x = elvm$model,  p=param, n = elvm$data$n, mu = elvm$mu, S = elvm$S)
  expect_equal(unname(Hred_num),unname(Hlava[indexRED,indexRED]), tolerance = 1e-6)
  
  Hred_I <- gaussianReduced_hessian.lvm(x = lvmRed, data=data, p=param[indexRED], type = "information")
  Hlava <- lava:::gaussian_hessian.lvm(x = elvm$model, data=data, p=param, n = elvm$data$n, mu = elvm$mu, S = elvm$S)
  # Hred_I-Hlava[indexRED,indexRED]

  return(invisible(list(Hred_E = Hred_E,
                        Hred_num = Hred_num,
                        Hred_I = Hred_I,
                        Hlava = Hlava)))
}



# gaussian_gradient.lvm <- function (x, p, data, S, mu, n, ...) 
# {
#   dots <- list(...)
#   dots$weight <- NULL
#   if (n > 2) 
#     data <- NULL
#   val <- -gaussian_score.lvm(x, p = p, S = S, mu = mu, n = n, 
#                              data = data, reindex = FALSE, ...)
#   if (!is.null(nrow(val))) {
#     val <- colSums(val)
#   }
#   val
# }


estimate.lvm <- function (x, data = parent.frame(), estimator = "gaussian", control = list(), 
          missing = FALSE, weight, weightname, weight2, id, fix, index = TRUE, 
          graph = FALSE, silent = lava.options()$silent, quick = FALSE, 
          method, param, cluster, p, ...) 
{
  if (length(exogenous(x) > 0)) {
    catx <- categorical2dummy(x, data)
    x <- catx$x
    data <- catx$data
  }
  cl <- match.call()
  if (!base::missing(param)) {
    oldparam <- lava.options()$param
    lava.options(param = param)
    on.exit(lava.options(param = oldparam))
  }
  if (!base::missing(method)) {
    control["method"] <- list(method)
  }
  optim <- list(iter.max = lava.options()$iter.max, trace = ifelse(lava.options()$debug, 
                                                                   3, 0), gamma = lava.options()$gamma, gamma2 = 1, ngamma = lava.options()$ngamma, 
                lambda = 0.05, abs.tol = 1e-09, epsilon = 1e-10, delta = 1e-10, 
                rel.tol = 1e-10, S.tol = 1e-05, stabil = FALSE, start = NULL, 
                constrain = lava.options()$constrain, method = NULL, 
                starterfun = "startvalues0", information = "E", meanstructure = TRUE, 
                sparse = FALSE, tol = lava.options()$tol)
  defopt <- lava.options()[]
  defopt <- defopt[intersect(names(defopt), names(optim))]
  optim[names(defopt)] <- defopt
  if (length(control) > 0) {
    optim[names(control)] <- control
  }
  if (is.environment(data)) {
    innames <- intersect(ls(envir = data), vars(x))
    data <- as.data.frame(lapply(innames, function(x) get(x, 
                                                          envir = data)))
    names(data) <- innames
  }
  if (!lava.options()$exogenous) 
    exogenous(x) <- NULL
  redvar <- intersect(intersect(parlabels(x), latent(x)), colnames(data))
  if (length(redvar) > 0) 
    warning(paste("Latent variable exists in dataset", redvar))
  xfix <- setdiff(colnames(data)[(colnames(data) %in% parlabels(x, 
                                                                exo = TRUE))], latent(x))
  if (base::missing(fix)) {
    fix <- ifelse(length(xfix) > 0, FALSE, TRUE)
  }
  Debug(list("start=", optim$start))
  if (!base::missing(cluster)) 
    id <- cluster
  if (!base::missing(weight)) {
    if (is.character(weight)) {
      weight <- data[, weight, drop = FALSE]
      if (!base::missing(weightname)) {
        colnames(weight) <- weightname
      }
      else {
        yvar <- index(x)$endogenous
        nw <- seq_len(min(length(yvar), ncol(weight)))
        colnames(weight)[nw] <- yvar[nw]
      }
    }
    weight <- cbind(weight)
  }
  else {
    weight <- NULL
  }
  if (!base::missing(weight2)) {
    if (is.character(weight2)) {
      weight2 <- data[, weight2]
    }
  }
  else {
    weight2 <- NULL
  }
  if (!base::missing(id)) {
    if (is.character(id)) {
      id <- data[, id]
    }
  }
  else {
    id <- NULL
  }
  Debug("procdata")
  dd <- procdata.lvm(x, data = data)
  S <- dd$S
  mu <- dd$mu
  n <- dd$n
  {
    var.missing <- setdiff(vars(x), colnames(S))
    if (length(var.missing) > 0) {
      new.lat <- setdiff(var.missing, latent(x))
      if (length(new.lat) > 0) 
        x <- latent(x, new.lat)
    }
  }
  myhooks <- gethook()
  for (f in myhooks) {
    res <- do.call(f, list(x = x, data = data, weight = weight, 
                           weight2 = weight2, estimator = estimator, optim = optim))
    if (!is.null(res$x)) 
      x <- res$x
    if (!is.null(res$data)) 
      data <- res$data
    if (!is.null(res$weight)) 
      weight <- res$weight
    if (!is.null(res$weight2)) 
      weight2 <- res$weight2
    if (!is.null(res$optim)) 
      optim <- res$optim
    if (!is.null(res$estimator)) 
      estimator <- res$estimator
    rm(res)
  }
  checkestimator <- function(x, ...) {
    ffname <- paste0(x, c("_objective", "_gradient"), ".lvm")
    exists(ffname[1]) || exists(ffname[2])
  }
  if (!checkestimator(estimator)) {
    estimator <- tolower(estimator)
    if (!checkestimator(estimator)) {
      estimator <- toupper(estimator)
    }
  }
  ObjectiveFun <- paste0(estimator, "_objective", ".lvm")
  GradFun <- paste0(estimator, "_gradient", ".lvm")
  if (!exists(ObjectiveFun) & !exists(GradFun)) 
    stop("Unknown estimator.")
  Method <- paste0(estimator, "_method", ".lvm")
  if (!exists(Method)) {
    Method <- "nlminb1"
  }
  else {
    Method <- get(Method)
  }
  NoOptim <- "method" %in% names(control) && is.null(control$method)
  if (is.null(optim$method) && !(NoOptim)) {
    optim$method <- if (missing) 
      "nlminb1"
    else Method
  }
  if (!quick & index) {
    x <- fixsome(x, measurement.fix = fix, S = S, mu = mu, 
                 n = n, debug = !silent)
    if (!silent) 
      message("Reindexing model...\n")
    if (length(xfix) > 0) {
      index(x) <- reindex(x, sparse = optim$sparse, zeroones = TRUE, 
                          deriv = TRUE)
    }
    else {
      x <- updatelvm(x, sparse = optim$sparse, zeroones = TRUE, 
                     deriv = TRUE, mean = TRUE)
    }
  }
  if (is.null(estimator) || estimator == FALSE) {
    return(x)
  }
  if (length(index(x)$endogenous) == 0) 
    stop("No observed outcome variables. Check variable names in model and data.")
  if (!optim$meanstructure) {
    mu <- NULL
  }
  nparall <- index(x)$npar + ifelse(optim$meanstructure, index(x)$npar.mean + 
                                      index(x)$npar.ex, 0)
  if (!missing(p)) {
    start <- p
    optim$start <- p
  }
  else {
    myparnames <- coef(x, mean = TRUE)
    paragree <- FALSE
    paragree.2 <- c()
    if (!is.null(optim$start)) {
      paragree <- myparnames %in% names(optim$start)
      paragree.2 <- names(optim$start) %in% myparnames
    }
    if (sum(paragree) >= length(myparnames)) 
      optim$start <- optim$start[which(paragree.2)]
    if (!(length(optim$start) == length(myparnames) & sum(paragree) == 
          0)) 
      if (is.null(optim$start) || sum(paragree) < length(myparnames)) {
        if (is.null(optim$starterfun) && lava.options()$param != 
            "relative") 
          optim$starterfun <- startvalues0
        start <- suppressWarnings(do.call(optim$starterfun, 
                                          list(x = x, S = S, mu = mu, debug = lava.options()$debug, 
                                               silent = silent, data = data, ...)))
        if (!is.null(x$expar) && length(start) < nparall) {
          ii <- which(index(x)$e1 == 1)
          start <- c(start, structure(unlist(x$expar[ii]), 
                                      names = names(x$expar)[ii]))
        }
        if (length(paragree.2) > 0) {
          start[which(paragree)] <- optim$start[which(paragree.2)]
        }
        optim$start <- start
      }
  }
  coefname <- coef(x, mean = optim$meanstructure, fix = FALSE)
  names(optim$start) <- coefname
  if (missing) {
    return(estimate.MAR(x = x, data = data, fix = fix, control = optim, 
                        debug = lava.options()$debug, silent = silent, estimator = estimator, 
                        weight = weight, weight2 = weight2, cluster = id, 
                        ...))
  }
  constr <- lapply(constrain(x), function(z) (attributes(z)$args))
  xconstrain <- intersect(unlist(constr), manifest(x))
  xconstrainM <- TRUE
  XconstrStdOpt <- TRUE
  if (length(xconstrain) > 0) {
    constrainM <- names(constr) %in% unlist(x$mean)
    for (i in seq_len(length(constr))) {
      if (!constrainM[i]) {
        if (constr[[i]] %in% xconstrain) {
          xconstrainM <- FALSE
          break
        }
      }
    }
    if (xconstrainM & ((is.null(control$method) || optim$method == 
                        "nlminb0") & (lava.options()$test & estimator == 
                                      "gaussian"))) {
      XconstrStdOpt <- FALSE
      optim$method <- "nlminb0"
      if (is.null(control$constrain)) 
        control$constrain <- TRUE
    }
  }
  lowmin <- -Inf
  lower <- rep(lowmin, length(optim$start))
  if (length(optim$constrain) == 1 & optim$constrain) 
    lower[variances(x) + index(x)$npar.mean] <- optim$tol
  if (any(optim$constrain)) {
    if (length(optim$constrain) != length(lower)) 
      constrained <- is.finite(lower)
    else constrained <- optim$constrain
    lower[] <- -Inf
    optim$constrain <- TRUE
    constrained <- which(constrained)
    nn <- names(optim$start)
    CS <- optim$start[constrained]
    CS[CS < 0] <- 0.01
    optim$start[constrained] <- log(CS)
    names(optim$start) <- nn
  }
  optim$start[is.nan(optim$start)] <- 0
  ObjectiveFun <- paste0(estimator, "_objective", ".lvm")
  GradFun <- paste0(estimator, "_gradient", ".lvm")
  if (!exists(ObjectiveFun) & !exists(GradFun)) 
    stop("Unknown estimator.")
  InformationFun <- paste0(estimator, "_hessian", ".lvm")
  mymodel <- x
  myclass <- "lvmfit"
  if (length(xfix) > 0 | (length(xconstrain) > 0 & XconstrStdOpt | 
                          !lava.options()$test)) {
    x0 <- x
    if (length(xfix) > 0) {
      myclass <- c("lvmfit.randomslope", myclass)
      nrow <- length(vars(x))
      xpos <- lapply(xfix, function(y) which(regfix(x)$labels == 
                                               y))
      colpos <- lapply(xpos, function(y) ceiling(y/nrow))
      rowpos <- lapply(xpos, function(y) (y - 1)%%nrow + 
                         1)
      myfix <- list(var = xfix, col = colpos, row = rowpos)
      x0 <- x
      for (i in seq_along(myfix$var)) for (j in seq_len(length(myfix$col[[i]]))) regfix(x0, 
                                                                                        from = vars(x0)[myfix$row[[i]][j]], to = vars(x0)[myfix$col[[i]][j]]) <- colMeans(data[, 
                                                                                                                                                                               myfix$var[[i]], drop = FALSE])
      x0 <- updatelvm(x0, zeroones = TRUE, deriv = TRUE)
      x <- x0
      yvars <- endogenous(x0)
      new.par.idx <- which(coef(mymodel, mean = TRUE, fix = FALSE) %in% 
                             coef(x0, mean = TRUE, fix = FALSE))
      if (length(optim$start) > length(new.par.idx)) 
        optim$start <- optim$start[new.par.idx]
      lower <- lower[new.par.idx]
      if (optim$constrain) {
        constrained <- match(constrained, new.par.idx)
      }
    }
    mydata <- as.matrix(data[, manifest(x0)])
    myObj <- function(pp) {
      if (optim$constrain) {
        pp[constrained] <- exp(pp[constrained])
      }
      myfun <- function(ii) {
        if (length(xfix) > 0) 
          for (i in seq_along(myfix$var)) {
            x0$fix[cbind(rowpos[[i]], colpos[[i]])] <- index(x0)$A[cbind(rowpos[[i]], 
                                                                         colpos[[i]])] <- data[ii, xfix[i]]
          }
        if (is.list(weight2)) {
          res <- do.call(ObjectiveFun, list(x = x0, p = pp, 
                                            data = mydata[ii, ], n = 1, weight = weight[ii, 
                                                                                        ], weight2 = weight2[ii, ]))
        }
        else {
          res <- do.call(ObjectiveFun, list(x = x0, p = pp, 
                                            data = mydata[ii, ], n = 1, weight = weight[ii, 
                                                                                        ], weight2 = weight2))
        }
        return(res)
      }
      sum(sapply(seq_len(nrow(data)), myfun))
    }
    myGrad <- function(pp) {
      if (optim$constrain) {
        pp[constrained] <- exp(pp[constrained])
      }
      myfun <- function(ii) {
        if (length(xfix) > 0) 
          for (i in seq_along(myfix$var)) {
            x0$fix[cbind(rowpos[[i]], colpos[[i]])] <- index(x0)$A[cbind(rowpos[[i]], 
                                                                         colpos[[i]])] <- data[ii, xfix[i]]
          }
        if (is.list(weight2)) {
          rr <- do.call(GradFun, list(x = x0, p = pp, 
                                      data = mydata[ii, , drop = FALSE], n = 1, 
                                      weight = weight[ii, ], weight2 = weight2))
        }
        else {
          rr <- do.call(GradFun, list(x = x0, p = pp, 
                                      data = mydata[ii, , drop = FALSE], n = 1, 
                                      weight = weight[ii, ], weight2 = weight2[ii, 
                                                                               ]))
        }
        return(rr)
      }
      ss <- rowSums(rbind(sapply(seq_len(nrow(data)), myfun)))
      if (optim$constrain) {
        ss[constrained] <- ss[constrained] * pp[constrained]
      }
      return(ss)
    }
    myInfo <- function(pp) {
      myfun <- function(ii) {
        if (length(xfix) > 0) 
          for (i in seq_along(myfix$var)) {
            x0$fix[cbind(rowpos[[i]], colpos[[i]])] <- index(x0)$A[cbind(rowpos[[i]], 
                                                                         colpos[[i]])] <- data[ii, xfix[i]]
          }
        if (is.list(weight2)) {
          res <- do.call(InformationFun, list(p = pp, 
                                              obj = myObj, x = x0, data = data[ii, ], n = 1, 
                                              weight = weight[ii, ], weight2 = weight2))
        }
        else {
          res <- do.call(InformationFun, list(p = pp, 
                                              obj = myObj, x = x0, data = data[ii, ], n = 1, 
                                              weight = weight[ii, ], weight2 = weight2[ii, 
                                                                                       ]))
        }
        return(res)
      }
      L <- lapply(seq_len(nrow(data)), function(x) myfun(x))
      val <- apply(array(unlist(L), dim = c(length(pp), 
                                            length(pp), nrow(data))), c(1, 2), sum)
      if (!is.null(attributes(L[[1]])$grad)) {
        attributes(val)$grad <- colSums(matrix(unlist(lapply(L, 
                                                             function(i) attributes(i)$grad)), ncol = length(pp), 
                                               byrow = TRUE))
      }
      return(val)
    }
  }
  else {
    xconstrain <- c()
    for (i in seq_len(length(constrain(x)))) {
      z <- constrain(x)[[i]]
      xx <- intersect(attributes(z)$args, manifest(x))
      if (length(xx) > 0) {
        warg <- setdiff(attributes(z)$args, xx)
        wargidx <- which(attributes(z)$args %in% warg)
        exoidx <- which(attributes(z)$args %in% xx)
        parname <- names(constrain(x))[i]
        y <- names(which(unlist(lapply(intercept(x), 
                                       function(x) x == parname))))
        el <- list(i, y, parname, xx, exoidx, warg, wargidx, 
                   z)
        names(el) <- c("idx", "endo", "parname", "exo", 
                       "exoidx", "warg", "wargidx", "func")
        xconstrain <- c(xconstrain, list(el))
      }
    }
    yconstrain <- unlist(lapply(xconstrain, function(x) x$endo))
    iconstrain <- unlist(lapply(xconstrain, function(x) x$idx))
    MkOffset <- function(pp, grad = FALSE) {
      if (length(xconstrain) > 0) {
        Mu <- matrix(0, nrow(data), length(vars(x)))
        colnames(Mu) <- vars(x)
        M <- modelVar(x, p = pp, data = data)
        M$parval <- c(M$parval, x$mean[unlist(lapply(x$mean, 
                                                     is.numeric))])
        for (i in seq_len(length(xconstrain))) {
          pp <- unlist(M$parval[xconstrain[[i]]$warg])
          myidx <- with(xconstrain[[i]], order(c(wargidx, 
                                                 exoidx)))
          mu <- with(xconstrain[[i]], apply(data[, exo, 
                                                 drop = FALSE], 1, function(x) func(unlist(c(pp, 
                                                                                             x))[myidx])))
          Mu[, xconstrain[[i]]$endo] <- mu
        }
        offsets <- Mu %*% t(M$IAi)[, endogenous(x)]
        return(offsets)
      }
      return(NULL)
    }
    myObj <- function(pp) {
      if (optim$constrain) {
        pp[constrained] <- exp(pp[constrained])
      }
      offset <- MkOffset(pp)
      mu0 <- mu
      S0 <- S
      x0 <- x
      if (!is.null(offset)) {
        x0$constrain[iconstrain] <- NULL
        data0 <- data[, manifest(x0)]
        data0[, endogenous(x)] <- data0[, endogenous(x)] - 
          offset
        pd <- procdata.lvm(x0, data = data0)
        S0 <- pd$S
        mu0 <- pd$mu
        x0$mean[yconstrain] <- 0
      }
      do.call(ObjectiveFun, list(x = x0, p = pp, data = data, 
                                 S = S0, mu = mu0, n = n, weight = weight, weight2 = weight2, 
                                 offset = offset))
    }
    myGrad <- function(pp) {
      if (optim$constrain) 
        pp[constrained] <- exp(pp[constrained])
      S <- do.call(GradFun, list(x = x, p = pp, data = data, 
                                 S = S, mu = mu, n = n, weight = weight, weight2 = weight2))
      if (optim$constrain) {
        S[constrained] <- S[constrained] * pp[constrained]
      }
      if (is.null(mu) & index(x)$npar.mean > 0) {
        return(S[-c(seq_len(index(x)$npar.mean))])
      }
      if (length(S) < length(pp)) 
        S <- c(S, rep(0, length(pp) - length(S)))
      return(S)
    }
    myInfo <- function(pp) {
      I <- do.call(InformationFun, list(p = pp, obj = myObj, 
                                        x = x, data = data, S = S, mu = mu, n = n, weight = weight, 
                                        weight2 = weight2, type = optim$information))
      if (is.null(mu) && index(x)$npar.mean > 0) {
        return(I[-seq_len(index(x)$npar.mean), -seq_len(index(x)$npar.mean)])
      }
      return(I)
    }
  }
  myHess <- function(pp) {
    p0 <- pp
    if (optim$constrain) 
      pp[constrained] <- exp(pp[constrained])
    I0 <- myInfo(pp)
    attributes(I0)$grad <- NULL
    D <- attributes(I0)$grad
    if (is.null(D)) {
      D <- myGrad(p0)
      attributes(I0)$grad <- D
    }
    if (optim$constrain) {
      I0[constrained, -constrained] <- apply(I0[constrained, 
                                                -constrained, drop = FALSE], 2, function(x) x * 
                                               pp[constrained])
      I0[-constrained, constrained] <- t(I0[constrained, 
                                            -constrained])
      if (sum(constrained) == 1) {
        I0[constrained, constrained] <- I0[constrained, 
                                           constrained] * outer(pp[constrained], pp[constrained]) - 
          D[constrained]
      }
      else {
        I0[constrained, constrained] <- I0[constrained, 
                                           constrained] * outer(pp[constrained], pp[constrained]) - 
          diag(D[constrained], ncol = length(constrained))
      }
    }
    return(I0)
  }
  if (is.null(tryCatch(get(InformationFun), error = function(x) NULL))) 
    myInfo <- myHess <- NULL
  if (is.null(tryCatch(get(GradFun), error = function(x) NULL))) 
    myGrad <- NULL
  if (!silent) 
    message("Optimizing objective function...")
  if (optim$trace > 0 & !silent) 
    message("\n")
  if ((is.data.frame(data) | is.matrix(data)) && nrow(data) == 
      0) 
    stop("No observations")
  if (!missing(p)) {
    opt <- list(estimate = p)
  }
  else {
    if (!is.null(optim$method)) {
      
      # [OB DEBUG]
      # myObj(optim$start)
      # setNames(myGrad(optim$start), names(optim$start))
      # myHess(optim$start)
      
      opt <- do.call(optim$method, list(start = optim$start, 
                                        objective = myObj, gradient = myGrad, hessian =  myHess,#hessian = myHess,   #
                                        lower = lower, control = optim, debug = debug))
      
      if (is.null(opt$estimate)) 
        opt$estimate <- opt$par
      if (optim$constrain) {
        opt$estimate[constrained] <- exp(opt$estimate[constrained])
      }
      if (XconstrStdOpt & !is.null(myGrad)) 
        opt$gradient <- as.vector(myGrad(opt$par))
      else {
        opt$gradient <- numDeriv::grad(myObj, opt$par)
      }
    }
    else {
      if (!NoOptim) {
        opt <- do.call(ObjectiveFun, list(x = x, data = data, 
                                          control = control, ...))
        opt$gradient <- rep(0, length(opt$estimate))
      }
      else {
        opt <- list(estimate = optim$start, gradient = rep(0, 
                                                           length(optim$start)))
      }
    }
  }
  if (!is.null(opt$convergence)) {
    if (opt$convergence != 0) 
      warning("Lack of convergence. Increase number of iteration or change starting values.")
  }
  else if (!is.null(opt$gradient) && mean(opt$gradient)^2 > 
           0.001) 
    warning("Lack of convergence. Increase number of iteration or change starting values.")
  if (quick) {
    return(opt$estimate)
  }
  pp <- rep(NA, length(coefname))
  names(pp) <- coefname
  if (!is.null(names(opt$estimate))) {
    pp[names(opt$estimate)] <- opt$estimate
    pp.idx <- na.omit(match(coefname, names(opt$estimate)))
  }
  else {
    pp[] <- opt$estimate
    pp.idx <- seq(length(pp))
  }
  suppressWarnings(mom <- tryCatch(modelVar(x, pp, data = data), 
                                   error = function(x) NULL))
  if (NoOptim) {
    asVar <- matrix(NA, ncol = length(pp), nrow = length(pp))
  }
  else {
    if (!silent) 
      message("\nCalculating asymptotic variance...\n")
    asVarFun <- paste0(estimator, "_variance", ".lvm")
    if (!exists(asVarFun)) {
      if (is.null(myInfo)) {
        if (!is.null(myGrad)) 
          myInfo <- function(pp) numDeriv::jacobian(myGrad, 
                                                    pp, method = lava.options()$Dmethod)
        else myInfo <- function(pp) numDeriv::hessian(myObj, 
                                                      pp)
      }
      I <- myInfo(opt$estimate)
      asVar <- tryCatch(Inverse(I), error = function(e) matrix(NA, 
                                                               length(opt$estimate), length(opt$estimate)))
    }
    else {
      asVar <- tryCatch(do.call(asVarFun, list(x = x, p = opt$estimate, 
                                               data = data, opt = opt)), error = function(e) matrix(NA, 
                                                                                                    length(opt$estimate), length(opt$estimate)))
    }
    if (any(is.na(asVar))) {
      warning("Problems with asymptotic variance matrix. Possibly non-singular information matrix!")
    }
    if (!is.null(attributes(asVar)$pseudo) && attributes(asVar)$pseudo) {
      warning("Near-singular covariance matrix, using pseudo-inverse!")
    }
    diag(asVar)[diag(asVar) == 0] <- NA
  }
  mycoef <- matrix(NA, nrow = nparall, ncol = 4)
  mycoef[pp.idx, 1] <- opt$estimate
  res <- list(model = x, call = cl, coef = mycoef, vcov = asVar, 
              mu = mu, S = S, model0 = mymodel, estimator = estimator, 
              opt = opt, expar = x$expar, data = list(model.frame = data, 
                                                      S = S, mu = mu, C = mom$C, v = mom$v, n = n, m = length(latent(x)), 
                                                      k = length(index(x)$manifest), weight2 = weight2), 
              weight = weight, weight2 = weight2, cluster = id, pp.idx = pp.idx, 
              graph = NULL, control = optim)
  class(res) <- myclass
  myhooks <- gethook("post.hooks")
  for (f in myhooks) {
    res0 <- do.call(f, list(x = res))
    if (!is.null(res0)) 
      res <- res0
  }
  if (graph) {
    res <- edgelabels(res, type = "est")
  }
  return(res)
}

NR <- function (start, objective, gradient, hessian, debug = FALSE, 
          control, ...) 
{ 
  control0 <- list(trace = 0, gamma = 1, lambda = 0, ngamma = 0, 
                   gamma2 = 0, backtrace = TRUE, iter.max = 200, tol = 1e-09, 
                   stabil = FALSE, epsilon = 1e-09)
  if (!missing(control)) {
    control0[names(control)] <- control
  }
  if (control0$trace > 0) 
    cat("\nIter=0;\t\n", "\tp=", paste0(formatC(start), collapse = " "), 
        "\n")
  gradFun = !missing(gradient)
  if (!gradFun & missing(hessian)) {
    hessian <- function(p) {
      ff <- objective(p)
      res <- attributes(ff)$hessian
      attributes(res)$grad <- as.vector(attributes(ff)$grad)
      return(res)
    }
  }
  oneiter <- function(p.orig, Dprev, return.mat = FALSE) {
    if (is.null(hessian)) {
      cat(".")
      I <- -numDeriv::jacobian(gradient, p.orig, method = lava.options()$Dmethod)
    }
    else {
      I <- -hessian(p.orig)
    }
    D <- attributes(I)$grad
    if (is.null(D)) {
      D <- gradient(p.orig)
    }
    print("****")
    print(D)
    print("")
    print(as.double(p.orig))
    print("****")
    if (return.mat) 
      return(list(D = D, I = I))
    if (control0$stabil) {
      if (control0$lambda != 0) {
        if (control0$lambda < 0) {
          sigma <- (t(D) %*% (D))[1]
        }
        else {
          sigma <- control0$lambda
        }
        sigma <- min(sigma, 10)
        I <- I + control0$gamma2 * sigma * diag(nrow = nrow(I))
      }
      else {
        sigma <- ((D) %*% t(D))
        I <- I + control0$gamma2 * (sigma)
      }
    }
    svdI <- svd(I)
    svdI$d0 <- numeric(length(svdI$d))
    svdI$d0[abs(svdI$d) > control0$epsilon] <- 1/svdI$d[abs(svdI$d) > 
                                                          control0$epsilon]
    iI <- with(svdI, (v) %*% diag(d0, nrow = length(d0)) %*% 
                 t(u))
    Delta = control0$gamma * iI %*% D
    Lambda <- 1
    if (control0$backtrace) {
      mD0 <- mean(Dprev^2)
      mD <- mean(D^2)
      p <- p.orig + Lambda * Delta
      while (mD >= mD0) {
        if (gradFun) {
          D = gradient(p)
        }
        else {
          DI <- oneiter(p, return.mat = TRUE)
          D = DI$D
        }
        mD = mean(D^2)
        if (is.nan(mD)) 
          mD = mD0
        Lambda <- Lambda/2
        if (Lambda < 1e-06) 
          break
        p <- p.orig + Lambda * Delta
      }
    }
    else {
      p <- p.orig + Lambda * Delta
    }
    return(list(p = p, D = D, iI = iI))
  }
  count <- count2 <- 0
  thetacur <- start
  gammacount <- 0
  Dprev <- rep(Inf, length(start))
  for (jj in seq_len(control0$iter.max)) {
    gammacount <- gammacount + 1
    count <- count + 1
    count2 <- count2 + 1
    oldpar <- thetacur
    newpar <- oneiter(thetacur, Dprev)
    Dprev <- newpar$D
    thetacur <- newpar$p
    if (!is.null(control0$ngamma) && control0$ngamma > 0) {
      if (control0$ngamma <= gammacount) {
        control0$gamma <- sqrt(control0$gamma)
        gammacount <- 0
      }
    }
    if (count2 == control0$trace) {
      cat("Iter=", count, ";\n\tD=", paste0(formatC(newpar$D), 
                                            collapse = " "), "\n")
      cat("\tp=", paste0(formatC(thetacur), collapse = " "), 
          "\n")
      count2 <- 0
    }
    if (mean(newpar$D^2) < control0$tol) 
      break
  }
  res <- list(par = as.vector(thetacur), iterations = count, 
              method = "NR", gradient = newpar$D, iH = newpar$iI)
  return(res)
}