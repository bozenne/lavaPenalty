APmatrix <- function(x){ # borrowed from coef.lvmfit
 
  names2.coef <- names(coef(x))
  if (is.null(x$control$meanstructure)){
    meanstructure <- TRUE
  } else {
    meanstructure <- x$control$meanstructure
  }
  npar <- index(x)$npar
  npar.mean <- index(x)$npar.mean*meanstructure
  npar.ex <- index(x)$npar.ex
  
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
      myorder.extra <- seq_len(index(x)$npar.ex)+length(myorder)
      myorder <- c(myorder,myorder.extra)
    }
  } else {
    myorder <- seq_len(npar+npar.mean)
    myorder.reg <- seq_len(npar)
    myorder.extra <- seq_len(index(x)$npar.ex)+length(myorder)
    myorder <- c(myorder,myorder.extra)
  }
  ## myorder <- seq_len(npar+npar.mean)
  ## myorder.reg <- seq_len(npar)
  ## myorder.extra <- seq_len(index(x)$npar.ex)+length(myorder)
  ## myorder <- c(myorder,myorder.extra)
  
  myparnames <- paste0("p",seq_len(npar+npar.ex))[myorder.reg]
  return(matrices.lvm(Model(x), myparnames))
  
}



#' @title Extract the type of the coefficient from a LVM
#' 
#' @param x a lvm model
#' 
#' @examples 
#' ## regression
#' m <- lvm(Y~X1+X2)
#' e <- estimate(m, sim(m, 1e2))
#' 
#' coefType(m)
#' coefType(e)
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
#' 
#' @export

#' @export
`coefType` <-
  function(x,...) UseMethod("coefType")

#' @export
coefType.lvm <- function(x, ...){ 
  
  ####
  names.coef <- coef(x)
  index.coef <- lava::index(x)
  
  type <- setNames(character(length = length(names.coef)), names.coef)
  type[index.coef$parBelongsTo$mean] <- "intercept"
  type[index.coef$parBelongsTo$reg] <- "regression"
  type[index.coef$parBelongsTo$cov] <- "covariance"
  type[index.coef$parBelongsTo$epar] <- "extra"
  
  #### variance
  type[names2.coef %in% diag(APmatrix(x)$P)] <- "variance"
  
  #### export
  return(type)
}

#' @export
coefType.lvmfit <- function(x, ...){ 
  
  ####
  type <- coefType.lvm(x$model0)
  
  #### export
  return(type)
}


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
  type <- coefType(x)
  names.var <- names(type)[type == "variance"]
  
  if(length(names.var)>0){
    return(grep(paste(names.var, collapse = "|"), coef(x), value = value))
  }else{
    return(NULL)
  }
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
  type <- coefType(x)
  
  validType <- "covariance"
  if(keep.var){
    validType <- c(validType, "variance")
  }
  
  names.cov <- names(type)[type %in% validType]
  if(length(names.cov)>0){
    return(grep(paste(names.cov, collapse = "|"), coef(x), value = value))
  }else{
    return(NULL)
  }
}


#' @title Extract the name or the position of the covariance coefficients
#' 
#' @param x a lvm model
#' @param value should the name of the coefficient be returned? Else return the coefficients
#' @param keep.var should the variance parameters be output?
#' 
#' @examples 
#' 
#' m <- lvm(Y~X1+X2)
#' categorical(m, K =3) <- "X1"
#' getIvar.lvm(m, "X1:0|1")
#' getIvar.lvm(m, c("X1:0|1", "X1:1|2"))
#' getIvar.lvm(m, "Y~X1")
#' getIvar.lvm(m, c("Y~X1","Y~X2"))
#' getIvar.lvm(m, c("Y~X2","Y~X1"))
#' 
getIvar.lvm <- function(x, link){
  
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
  
  #### deal with categorical variables
  if(length(index.Ncat)>0){
    link.Ncat <- link[index.Ncat]
    
    index(x)
    
    name.link <- names(coef(x))[match(link.Ncat,coef(x))]
    A <- APmatrix(x)$A
    colnames.A <- colnames(A)
    rownames.A <- rownames(A)
    
    M.Ncat <- sapply(name.link, function(l){
      position <- which(A == l, arr.ind = TRUE)
      return(c(colnames.A[position[2]],rownames.A[position[1]]))
    })
    colnames(M.Ncat) <- link
    rownames(M.Ncat) <- c("endogenous","exogenous")
  }else{
    M.Ncat <- NULL
  }
  
  #### deal with categorical variables
  if(length(index.cat)>0){
    
  link.cat <- link[index.cat]
  exo.link <- sapply(link.cat, function(l){
    test <- unlist(lapply(x$attributes$ordinalparname, function(vec.coef){l %in% vec.coef}))
    return(names(x$attributes$ordinalparname)[test])
  })
  
  if(any(rowSums(x$M[exo.link,,drop = FALSE])>1)){
    stop("Current implementation cannot deal with a categorical variable used for several outcomes \n")
  }
  
  endo.link <- apply(x$M[exo.link,,drop = FALSE], 1, function(row){ names(which(row==1)) })
  
  M.cat <- rbind(endo.link, exo.link)
  colnames(M.cat) <- link
  rownames(M.cat) <- c("endogenous","exogenous")
  }else{
    M.cat <- NULL
  }
  
  ## export
  return(cbind(M.Ncat,M.cat))
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