# NOTE: the user should not be able to remove the intercept from a model that have been reduced
`procdata` <- function(x,...) UseMethod("procdata") # to use procdata.lvm.reduced instead of procdata.lvm
deriv.lvm <- lava:::deriv.lvm # to use deriv.lvm without using :::

regression.lvm <- function(object = lvm(), to, from, fn = NA, silent = lava.options()$silent, 
                           additive = TRUE, reduce = FALSE, y, x, value, ...){
  
  if (!missing(y)) {
    if (inherits(y, "formula")) 
      y <- all.vars(y)
    to <- y
  }
  if (!missing(x)) {
    if (inherits(x, "formula")) 
      x <- all.vars(x)
    from <- x
  }
  if (!additive) {
    if (!inherits(to, "formula")) 
      to <- toformula(to, from)
    x <- attributes(getoutcome(to))$x
    K <- object$attributes$nordinal[x]
    if (is.null(K) || is.na(K)) {
      K <- list(...)$K
      if (is.null(K)) 
        stop("Supply number of categories, K (or use method 'categorical' before calling 'regression').")
      object <- categorical(object, x, ...)
    }
    dots <- list(...)
    dots$K <- K
    dots$x <- object
    dots$formula <- to
    dots$regr.only <- TRUE
    object <- do.call("categorical", dots)
    return(object)
  }
  if (missing(to)) {
    return(regfix(object))
  }
  if (inherits(to, "formula")) {
    if (!missing(value)) {
      regression(object, to, silent = silent, ...) <- value
    }
    else {
      regression(object, silent = silent, ...) <- to
    }
    object$parpos <- NULL
    return(object)
  }
  if (is.list(to)) {
    for (t in to) regression(object, silent = silent, ...) <- t
    object$parpos <- NULL
    return(object)
  }
  if(reduce){
    allCoef <- coef(regression(object, to = to, from = from, silent = silent, reduce = FALSE))
    
    if("lp" %in% names(object) == FALSE){object$lp <- list()}
    for(iterR in to){ ### I don't know where to get the constrains
      regression(object, to = iterR, from = paste0("LP",iterR), silent = silent) <- 1
      namesCoef <- paste(iterR, from, sep = "~")
      if(iterR %in% names(object$lp)){ 
        object$lp[[iterR]]$link <- c(object$lp[[iterR]]$link, namesCoef)
        # object$lp[[iterR]]$indexCoef <- match(object$lp[[iterR]]$indexCoef, coef(object))
        object$lp[[iterR]]$x <- c(object$lp[[iterR]]$x, from)
        object$lp[[iterR]]$con <- c(object$lp[[iterR]]$con, rep(NA,length(from)))
      }else{
        object$lp[[iterR]]$link <- namesCoef
        # object$lp[[iterR]]$indexCoef <- match(namesCoef, coef(object))
        object$lp[[iterR]]$x <- from
        object$lp[[iterR]]$con <- rep(NA,length(from))
        object$lp[[iterR]]$name <- paste0("LP",iterR)
      }
    }
    parameter(object) <- setdiff(allCoef,coef(object))
    class(object) <- append("lvm.reduced", class(object))
    return(object)
  }
  
  sx <- strsplit(from, "@")
  xx <- sapply(sx, FUN = function(i) i[1])
  ps <- sapply(sx, FUN = function(i) i[2])
  sx <- strsplit(xx, "$", fixed = TRUE)
  xs <- sapply(sx, FUN = function(i) i[1])
  fix <- as.numeric(sapply(sx, FUN = function(i) i[2]))
  allv <- index(object)$vars
  object <- addvar(object, c(to, xs), silent = silent, reindex = FALSE)
  for (i in to) for (j in xs) {
    object$M[j, i] <- 1
    if (!is.na(fn)) 
      functional(object, j, i) <- fn
  }
  if (lava.options()$exogenous) {
    newexo <- setdiff(xs, c(to, allv))
    exo <- exogenous(object)
    if (length(newexo) > 0) 
      exo <- unique(c(exo, newexo))
    exogenous(object) <- setdiff(exo, to)
  }
  if (lava.options()$debug) {
    print(object$fix)
  }
  object$fix[xs, to] <- fix
  object$par[xs, to] <- ps
  object$parpos <- NULL
  index(object) <- reindex(object)
  
  
  return(object)
}

