#' @title Extract the specific coefficient names or positions in a LVM
#' @name extractCoef
#' 
#' @param x a lvm model or a fitted lvm model 
#' @param value should the name of the coefficient be returned? Else return the coefficients
#' @param keep.var should the variance parameters be output?
#' 
#' @examples 
#' ## regression
#' m <- lvm(Y~X1+X2)
#' e <- estimate(m, sim(m, 1e2))
#' 
#' coefType(m)
#' coefType(e)
#' coefType(e, level = -1)
#'
#' coefCov(m)
#' coefCov(m, value = TRUE)
#'
#' coefCov(m, keep.var = TRUE)
#' coefCov(m, value = TRUE, keep.var = TRUE)
#'
#' coefIntercept(m)
#' coefIntercept(m, value = TRUE)
#'
#' coefReg(m)
#' coefReg(m, value = TRUE)
#' 
#' ## LVM
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x
#' latent(m) <- ~u
#' covariance(m) <- y1~y2
#' 
#' e <- estimate(m, sim(m, 1e2))
#' 
#' coefType(m)
#' coefType(e)
#' coefType(e, level = -1)
#'
#' coefCov(m)
#' coefCov(m, value = TRUE)#' 
#'
#' coefCov(m, keep.var = TRUE)
#' coefCov(m, value = TRUE, keep.var = TRUE)
#' 
#' coefExtra(m)
#'
#' categorical(m, K = 3) <- "X1"
#' coefExtra(m)
#' coefExtra(m, value = TRUE)
#'
#' coefIntercept(m)
#' coefIntercept(m, value = TRUE)
#' coefIntercept(e)
#'
#' coefReg(e, value = TRUE)
#' coefReg(e, level = -1, value = TRUE)
#' 


# {{{ coefType
#' @rdname extractCoef
#' @export
`coefType` <-
  function(x,...) UseMethod("coefType")

#' @rdname extractCoef
#' @export
coefType.lvm <- function(x, ...){ 
  
  ####
  names.coef <- coef(x, ...)
  index.coef <- lava::index(x)
  
  type <- setNames(character(length = length(names.coef)), names.coef)
  type[index.coef$parBelongsTo$mean] <- "intercept"
  type[index.coef$parBelongsTo$reg] <- "regression"
  type[index.coef$parBelongsTo$cov] <- "covariance"
  type[index.coef$parBelongsTo$epar] <- "extra"
  
  #### variance
  type[names(names.coef) %in% diag(APmatrix(x)$P)] <- "variance"
  
  #### export
  return(type)
}

#' @rdname extractCoef
#' @export
coefType.lvmfit <- function(x, level = 9, ...){ 
  
    #### all Coef
    extractCoef <- coef(x, level = level, ...)
    if(is.matrix(extractCoef)){
        names.allCoef <- rownames(extractCoef)
    }else{
        names.allCoef <- names(extractCoef)
    }

    res <- setNames(rep(NA, length = length(names.allCoef)), names.allCoef)
    attribute <- setNames(rep(TRUE, length = length(names.allCoef)), names.allCoef)

    #### with coef ref
    type <- coefType(x$model0, ...)        
    res[names(type)] <- type
    attribute[names(type)] <- FALSE

    ## try to recover coef ref type
    originalX <- try(eval(x$call$x))
    if("lvm" %in% class(originalX)){
        reftype <- coefType(originalX, ...)
        res[attribute==TRUE] <- reftype[names(attribute)[attribute==TRUE]]
    }   

    attr(res,"reference") <- attribute
        
    #### export
    return(res)
}
# }}}

# {{{ coefCov
#' @rdname extractCoef
#' @export
`coefCov` <-
  function(x,...) UseMethod("coefCov")

#' @rdname extractCoef
#' @export
coefCov.lvm <- function(x, value = FALSE, keep.var = FALSE, ...){

    res <- retainType(type = coefType(x, ...),
                      validType = c("covariance", if(keep.var){"variance"}else{NULL}),
                      value = value)

    return(res)
}

#' @rdname extractCoef
#' @export
coefCov.lvmfit <- coefCov.lvm

# }}}
# {{{ coefExtra
#' @rdname extractCoef
#' @export
`coefExtra` <-
  function(x,...) UseMethod("coefExtra")

#' @rdname extractCoef
#' @export
coefExtra.lvm <- function(x, value = FALSE, ...){ 

    res <- retainType(type = coefType(x, ...),
                      validType = "extra",
                      value = value) 
    
    return(res)    
}

#' @rdname extractCoef
#' @export
coefExtra.lvmfit <- coefExtra.lvm
# }}}
# {{{ coefIntercept
#' @rdname extractCoef
#' @export
`coefIntercept` <-
  function(x,...) UseMethod("coefIntercept")

#' @rdname extractCoef
#' @export
coefIntercept.lvm <- function(x, value = FALSE, ...){ 

    res <- retainType(type = coefType(x, ...),
                      validType = "intercept",
                      value = value)

    return(res)
}

#' @rdname extractCoef
#' @export
coefIntercept.lvmfit <- coefIntercept.lvm
# }}}
# {{{ coefRef
#' @rdname extractCoef
#' @export
`coefRef` <-
  function(x,...) UseMethod("coefRef")

#' @rdname extractCoef
#' @export
coefRef.lvmfit <- function(x, value = FALSE, ...){
    
    res <- retainType(type = attr(coefType(x, ...), "reference"),
                      validType = TRUE,
                      value = value)

    return(res)    
}

# }}}

# {{{ coefReg
#' @rdname extractCoef
#' @export
`coefReg` <-
  function(x,...) UseMethod("coefReg")

#' @rdname extractCoef
#' @export
coefReg.lvm <- function(x, value = FALSE, ...){
    
     res <- retainType(type = coefType(x, ...),
                      validType = "regression",
                      value = value)

     return(res)
}

#' @rdname extractCoef
#' @export
coefReg.lvmfit <- coefReg.lvm
# }}}

# {{{ coefVar
#' @rdname extractCoef
#' @export
`coefVar` <-
  function(x,...) UseMethod("coefVar")

#' @rdname extractCoef
#' @export
coefVar.lvm <- function(x, value = FALSE, ...){ 
    res <- retainType(type = coefType(x, ...),
                      validType = "variance",
                      value = value)

    return(res)
}

#' @rdname extractCoef
#' @export
coefVar.lvmfit <- coefVar.lvm
# }}}

# {{{ getIvar
#' @title Extract the variables related to each parameter
#' 
#' @param x a lvm model
#' @param link the links to be analysed. If NULL, all the coefficients from the lvm model are used instead.
#'
#' @param format the type of object to output. Can be \code{"data.frame"} or \code{"matrix"}.
#' 
#' @examples 
#' 
#' m <- lvm(Y~X1+X2)
#' categorical(m, K =3) <- "X1"
#' getIvar.lvm(m)
#' getIvar.lvm(m, "X1:0|1")
#' getIvar.lvm(m, c("X1:0|1", "X1:1|2"))
#' getIvar.lvm(m, "Y~X1")
#' getIvar.lvm(m, c("Y~X1","Y~X2"))
#' getIvar.lvm(m, c("Y~X2","Y~X1"))
#'
#' d <- sim(m, 1e2)
#' d$X1 <- factor(d$X1)
#' m2 <- lava_categorical2dummy(m, data = d)
#' getIvar.lvm(m2$x)
#' 
getIvar.lvm <- function(x, link = NULL, format = "data.frame"){

    if(format %in% c("data.frame","matrix") == FALSE){
        stop(paste0(format," is an invalid format. Can only be \"data.frame\" or \"matrix\" \n"))
    }

    #### valid links
    if(!is.null(link)){
        if(any(link %in% coef(x) == FALSE)){
            warning("unknown link(s): ", paste(link[link %in% coef(x) == FALSE], collpase = " "), "\n",
                    "possible link(s): ", paste(coef(x)[coef(x) %in% link == FALSE], collpase = " "), "\n")
            return(NULL)
        }
    }else{
        link <- coef(x)
    }

    index.cat <- which(link %in% unlist(x$attributes$ordinalparname))
    index.Ncat <- setdiff(1:length(link), index.cat)

  #### deal with continuous variables
  if(length(index.Ncat)>0){
    link.Ncat <- link[index.Ncat]
    
    name.link <- names(coef(x))[match(link.Ncat,coef(x))]
    A <- APmatrix(x)$A
    colnames.A <- colnames(A)
    rownames.A <- rownames(A)

    M.Ncat <- sapply(name.link, function(l){
      position <- which(A == l, arr.ind = TRUE)
      return(c(colnames.A[position[2]],rownames.A[position[1]]))
    })
    M.Ncat <- cbind(t(M.Ncat), "continous")
    
    rownames(M.Ncat) <- link.Ncat
    colnames(M.Ncat) <- c("endogenous","exogenous","type")
  }else{
    M.Ncat <- NULL
  }
  
  #### deal with categorical variables
  M.cat <- NULL
  if(length(index.cat)>0){
    
    link.cat <- link[index.cat]
    exo.link <- sapply(link.cat, function(l){
      test <- unlist(lapply(x$attributes$ordinalparname, function(vec.coef){l %in% vec.coef}))
      return(names(x$attributes$ordinalparname)[test])
    })
    
    # remove parameters corresponding to non existing variables (i.e. X1 when only X1B and X1C exists)
    exo.link <- exo.link[exo.link %in% rownames(x$M)]
    
    if(length(exo.link)>0){
      
      if(any(rowSums(x$M[exo.link,,drop = FALSE])>1)){
        stop("Current implementation cannot deal with a categorical variable used for several outcomes \n")
      }
      
      endo.link <- apply(x$M[exo.link,,drop = FALSE], 1, function(row){ names(which(row==1)) })
      M.cat <- cbind(endo.link, exo.link, "categorical")
      rownames(M.cat) <- link[index.cat]
      colnames(M.cat) <- c("endogenous","exogenous","type")
    }
    
  }
  
    ## export
    res <- rbind(M.Ncat, M.cat)
    if(format=="data.frame"){
        res <- as.data.frame(res, stringsAsFactors = FALSE)
    }    
    return(res)
}
# }}}

# {{{ loadings
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
#' @export
`loadings` <- function(object, ...) UseMethod("loadings")

#' @rdname loadings
#' @export
loadings.lvmfit <- function(x, level = 9, col = c("Estimate","Std. Error", "Z-value","P-value")){

    possibleCoef <- coefReg(x, level = level, value = TRUE)
    # check that the first variable is endogeneous and the second a latent variable
    test <- sapply(possibleCoef, function(var){
        vars <- initVar_link(var)
        test <- (vars[[1]] %in% endogenous(x))*(vars[[2]] %in% latent(x))
        return(rbind(test,vars[[1]],vars[[2]]))
    })

    tempo <- test[,test[1,] == "1"]
    names.loadings <- colnames(tempo)[order(tempo[3,])]
    
    # extraction 
    if(is.null(col)){
        return(names.loadings)
    }else{
        loadings <- summary(x)$coef[names.loadings,]
        match.arg(col, choices = c("Estimate","Std. Error", "Z-value","P-value"))
        return(loadings[,col])
    }
}

#' @rdname loadings
#' @export
loadings.lvm.missing <- loadings.lvmfit
# }}}

# {{{ Associated functions

# {{{ AP matrix (need for coefType)
APmatrix <- function(x){ # borrowed from coef.lvmfit
 
  names2.coef <- names(coef(x))
  if (is.null(x$control$meanstructure)){
    meanstructure <- TRUE
  } else {
    meanstructure <- x$control$meanstructure
  }
  npar <- lava::index(x)$npar
  npar.mean <- lava::index(x)$npar.mean*meanstructure
  npar.ex <- lava::index(x)$npar.ex
  
  if (inherits(x,"lvm.missing")) {
    if (length(x$cc)==0) {## No complete cases
      coefs <- coef(x$estimate)
      c1 <- coef(Model(x),mean=TRUE,fix=FALSE)
      c1. <- coef(Model(x),mean=FALSE,fix=FALSE)
      myorder <- match(c1,names(coefs))
      myorder.reg <- order(na.omit(match(names(coefs),c1.)))
      myorder.extra <- c()
      ##mp <-effect modelPar(x,seq_len(npar+npar.mean+npar.ex))
      ## mp <- modelPar(x,seq_len(npar+npar.mean+npar.ex))
      ## myorder <- c(mp$meanpar,mp$p)
      ## myorder.reg <- seq_len(length(mp$p))
      ## myorder.extra <- mp$p2
    } else {
      myorder <- na.omit(modelPar(x$multigroup,seq_len(npar+npar.mean))$p[[x$cc]])
      myorder.reg <- na.omit(modelPar(x$multigroup,seq_len(npar))$p[[x$cc]])
      myorder.extra <- seq_len(lava::index(x)$npar.ex)+length(myorder)
      myorder <- c(myorder,myorder.extra)
    }
  } else {
    myorder <- seq_len(npar+npar.mean)
    myorder.reg <- seq_len(npar)
    myorder.extra <- seq_len(lava::index(x)$npar.ex)+length(myorder)
    myorder <- c(myorder,myorder.extra)
  }
  ## myorder <- seq_len(npar+npar.mean)
  ## myorder.reg <- seq_len(npar)
  ## myorder.extra <- seq_len(lava::index(x)$npar.ex)+length(myorder)
  ## myorder <- c(myorder,myorder.extra)
  
  myparnames <- paste0("p",seq_len(npar+npar.ex))[myorder.reg]
  return(lava_matrices.lvm(Model(x), myparnames))
  
}
# }}}


# {{{ retainType (need for coefCov/Latent/Ref/Ref)
retainType <- function(type, validType, value){
  index.var <- which(type %in% validType)
  
  if(length(index.var)>0){
      if(value){
          return(names(type)[index.var])
      }else{
          return(index.var)
      }
  }else{
      return(NULL)
  }
}

# }}}

# }}}

# {{{ checkLatent
#' @title Extract the name or the position of the variance coefficients
#' @name checkVar
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
#' d <- sim(m,1e2)
#' checkVar(m, d)
#' 
#' checkVar(m, d[,-4])
#' @export
`checkVar` <-
  function(x,data,...) UseMethod("checkVar")

#' @rdname checkVar
#' @export
checkVar.lvm <- function(x, data){ 
    vars <- vars(x)
    latent <- latent(x)
    missingVars <- vars[vars %in% names(data) == FALSE]

    if(!identical(sort(latent),sort(missingVars))){
        cat("wrong specification of the latent variables \n",
            "latent variables according to the LVM: ",paste(latent, collapse = " "),"\n",
            "missing variables in data: ",paste(missingVars, collapse = " "),"\n")
        return(invisible(FALSE))
    }else{
        return(invisible(TRUE))
    }
}
# }}}
