#### intialisation ####

#' @title Initialize penalty terms
#' @name initPenalty
#' @description Initialise lasso, ridge, group lasso, and nuclear norm penalties
#' 
#' @details 
#' initPenaltyL12 contains the specification of the lasso, ridge, group lasso penalties
#' initPenaltNuclear contains the specification of the nuclear norm penalty

#' @rdname initPenalty
initPenaltyL12 <- function(){
  
  penalty <- list(lambda1 = 0, 
                  lambda2 = 0,
                  adaptive = FALSE,
                  proxOperator = NULL,
                  objectivePenalty = NULL,
                  Vlasso = NULL,
                  Vridge = NULL,
                  Vgroup = NULL)
  
  class(penalty) <- "penaltyL12"
  return(penalty)
  
}

#' @rdname initPenalty
initPenaltNuclear <- function(){
  
  penaltyNuclear <- list(link = NULL,
                         lambdaN  = 0, 
                         name.Y = NULL,
                         name.X = NULL,
                         nrow = NULL,
                         ncol = NULL)
  
  class(penaltyNuclear) <- "penaltyNuclear"
  return(penaltyNuclear)
  
}

#### convertion ####

#' @title Conversion to a penalized latent variable model
#' 
#' @param x \code{lvm}-object
#' 
lvm2plvm <- function(x){
  
  if("plvm" %in% class(x) == FALSE){
    x$penalty <- initPenaltyL12()
    x$penaltyNuclear <- initPenaltNuclear()
    class(x) <- append("plvm", class(x))
  }else{
    warning("x is already a penalized latent variable model \n")
  }
  
  return(x)
}

#### specification ####

#' @title Penalize a latent variable model
#' @description Add a penalty term to a latent variable model
#' @name penalize
#' 
#' @param x \code{lvm}-object
#' @param value the name of the link to be penalized 
#' @param group the groups defining the group lasso penalty
#' @parma coords the (spatial) position of the links to penalize
#' @param V the matrix defining lasso penalties
#' @param add should value be added to the existing penalty term ? Otherwise it will overwrite it.
#' @param reduce should for each regression the penalised link be aggregated into a linear predictor.
#' @param lambda1 lasso penalization parameter
#' @param lambda2 ridge penalization parameter
#' @param lambdaN nuclear norm penalization parameter
#' @param adaptive should an adaptive lasso be used?
#' @param intercept should all intercept be penalized. Disregarded if value is specified.
#' @param regression should all regression parameters be penalized. Disregarded if value is specified.
#' @param variance should all covariance links be penalized. Disregarded if value is specified.
#' @param latent If FALSE, no link related to the latent variable will be penalized. Disregarded if value is specified.
#'
#' @details 
#' By default categorical variables are penalised using a group lasso penalty.
#' penalize functions can be used to add lasso or/and ridge penalty terms to the model.
#' penaltyNuclear functions can be used to add a nuclear penalty to the model.
#' 
#' @examples 
#' 
#' #### lasso penalty ####
#' set.seed(10)
#' n <- 500
#' formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
#' lvm.modelSim <- lvm()
#' regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
#' distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
#' df.data <- sim(lvm.modelSim,n)
#' 
#' lvm.model <- lvm(formula.lvm)
#' plvm.model <- penalize(lvm.model)
#'
#' #### group penalty ####
#' m <- lvm(Y1 ~ X1+X2+X3+X4)
#' categorical(m, labels = c("A","B","C")) <- "X1"
#' categorical(m, labels = c("A","B","C")) <- "X2"
#' 
#' pm <- penalize(m, value = c("Y1~X1","Y1~X3","Y1~X4"))
#' pm$penalty
#' pm1 <- penalize(m, value = c("Y1~X1"))
#' pm1 <- penalize(pm1, value = c("Y1~X2","Y1~X3"))
#' 
#'
#'
#'m <- lvm(list(Y1 ~ X1+X2+X3+X4+eta, Y2 ~ X1+X2+X3+eta, Y3 ~ eta))
#'categorical(m, labels = c("A","B","C")) <- "X1"
#'latent(m) <- ~eta
#'
#'penalize(m)
#'


#'
#' @export
`penalize` <-
  function(x,...) UseMethod("penalize")

#' @name penalize
#' @export
"penalize<-" <- function (x, ..., value) {
  UseMethod("penalize<-", x)
}

#' @name penalize
#' @export
`penalize.lvm` <- function(x, value = NULL, ...){
  
  penalize(x, ...) <- value
  
  return(x)
}

#' @name penalize
#' @export
`penalize<-.lvm` <- function(x, group, V, add = TRUE, reduce = FALSE,
                             lambda1, lambda2, adaptive,
                             intercept = FALSE, regression = TRUE, variance = FALSE, covariance = FALSE, latent = FALSE,
                             value){
  
  ## convert to plvm
  if("plvm" %in% class(x) == FALSE){
    x <- lvm2plvm(x)
  }
  
  #### find coefficients from value
  if(!is.null(value)){
    
    if("formula" %in% class(value)){
      value <- formula2character(value)
    }else if(is.list(value)){
      value <- unlist(lapply(value, function(v){
        if("formula" %in% class(v)){formula2character(v)}else{v}
      }))
    }
    
    if(any(value %in% coef(x) == FALSE)){
      stop("penalty<-.lvm: coefficients to be penalized do not match those of the model\n",
           "unknown coefficients: ",paste(value[value %in% coef(x) == FALSE], collapse = " "),"\n",
           "available coefficients: ",paste(coef(x)[coef(x) %in% value == FALSE], collapse = " "),"\n")
    }
    
  } else if(is.null(penalty(x, type = "link"))){
    
    index.penaltyCoef <- NULL
    if(intercept == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, x$index$parBelongsTo$mean)  
    }
    if(regression == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, x$index$parBelongsTo$reg) 
    }
    if(variance == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, coefVar(x))
    }
    
    if(covariance == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, coefCov(x))
    }
    
    ## no penalization on parameters related to the latent variables
    if(length(x$latent) && latent == FALSE){
      test.latent <- sapply(coef(x)[index.penaltyCoef], function(pen){initVar_link(pen)$var1 %in% latent(x)})
      index.penaltyCoef <- setdiff(index.penaltyCoef, index.penaltyCoef[which(test.latent)])
    }
    
    value <- coef(x)[index.penaltyCoef]
  } 
  
  #### update the V matrices
  if(!missing(V)){
    if("Matrix" %in% is(V)){
      stop("V must herit from the class Matrix \n")
    }
    penalty(x, type = "V") <- V
  }else{
    resInit <- initGroup.lvm(x, links = value, group = group)
    
    test.groupPenalty <- !is.null(resInit$Mcat)
    
    if(test.groupPenalty){ #### Vgroup
      Vgroup.old <- penalty(x, type = "Vgroup")
      
      if(add && !is.null(Vgroup.old)){
        
        newV <- initVcoef.lvm(x, link = colnames(resInit$Mcat), group = resInit$Mcat["group",])
        newV@x <- newV@x + max(Vgroup.old@x)
        penalty(x, type = "Vgroup") <- cbind(Vgroup.old,newV)
        
      }else{
        penalty(x, type = "Vgroup") <- initVcoef.lvm(x, link = colnames(resInit$Mcat), group = resInit$Mcat["group",])
      }
    }
    
    if(NCOL(resInit$M)>0){ #### lasso
      penalty(x, type = "Vlasso") <- initVcoef.lvm(x, link = colnames(resInit$M), group = rep(1,NCOL(resInit$M)))  
    }
    
    if(NCOL(resInit$M)>0 || NCOL(resInit$Mcat)>0){ #### ridge
      penalty(x, type = "Vridge") <- initVcoef.lvm(x, link = c(colnames(resInit$M),colnames(resInit$Mcat)), group = rep(1,NCOL(resInit$M)+NCOL(resInit$Mcat)))
    }
    
  }
  
  #### penalization parameters
  if(!missing(lambda1)){
    penalty(x, type = "lambda1") <- as.numeric(lambda1)
  }
  
  if(!missing(lambda2)){
    penalty(x, type = "lambda2") <- as.numeric(lambda2)
  }
  
  if(!missing(adaptive)){
    penalty(x, type = "adaptive") <- as.numeric(adaptive)
  }
  
  #### reduce 
  if(reduce){
    x <- reduce(x)
  }
  
  
  #### export
  return(x)
}

#' @name penalize
#' @export
`penalizeNuclear` <-
  function(x,...) UseMethod("penalizeNuclear")

#' @name penalize
#' @export
"penalizeNuclear<-" <- function (x, ..., value) {
  UseMethod("penalizeNuclear<-", x)
}

#' @name penalize
#' @export
`penalizeNuclear.lvm` <- function(x, value = NULL, ...){
  
  penalizeNuclear(x, ...) <- value
  
  return(x)
}

#' @name penalize
#' @export
`penalizeNuclear<-.plvm` <- function(x, coords, lambdaN = NULL, ..., value){
  
  ## convert to plvm
  if("plvm" %in% class(x) == FALSE){
    x <- lvm2plvm(x)
  }
  
  ## lambda
  if(!is.null(lambdaN)){
    penalty(x, type = "lambdaN", nuclear = TRUE) <- lambdaN
  }
  
  ## formula
  name.Y <- all.vars(value)[1]
  name.X <- all.vars(value)[-1]
  
  if(name.Y %in% endogenous(x) == FALSE){
    stop("penaltyNuclear: the dependent variable in formula must already be in the model \n")
  }
  if(any(name.X %in% vars(x) == TRUE)){
    stop("penaltyNuclear: the independent variable in formula must not be in the model \n",
         "existing variables: ",paste(name.X[name.X %in% vars(x) == TRUE],collapse = " "),"\n")
  }
  coords.factor <- apply(coords, 2, function(x){
    as.numeric(as.factor(x))
  })
  ncol <- max(coords.factor)
  nrow <- max(coords.factor)
  
  if(ncol*nrow != NROW(coords)){
    stop("penaltyNuclear: coords does not corresponds to a full matrix \n")
  }
  test.row <- sweep(matrix(coords.factor[,1], nrow = nrow, ncol = ncol),
                    MARGIN = 1, STATS = 1:nrow, FUN = "-")
  if(any(test.row!=0) ){
    stop("penaltyNuclear: coordinates must be ordered \n")
  }
  test.col <- sweep(matrix(coords.factor[,2], nrow = nrow, ncol = ncol),
                    MARGIN = 2, STATS = 1:ncol, FUN = "-")
  if(any(test.col!=0) ){
    stop("penaltyNuclear: coordinates must be ordered \n")
  }
  
  #### update penalty
  penalty(x, type = "link", nuclear = TRUE) <- paste0(name.Y,"~",name.X)
  penalty(x, type = "name.Y", nuclear = TRUE) <- name.Y
  penalty(x, type = "name.X", nuclear = TRUE) <- name.X
  penalty(x, type = "ncol", nuclear = TRUE) <- ncol
  penalty(x, type = "nrow", nuclear = TRUE) <- nrow
  
  #### update object
  x <- regression(x, from = name.X, to = name.Y, reduce = paste0("Image_",LCSseq(name.X),"_")) 
  
  #### export
  return(x)
}


#### init ####

#' @title Reshape information for building the V matrices
#' @description Return matrices containg the endogenous, exogenous and group index corresponding to each penlized link. Take care of categorical variables.
#' 
#' @param x a \code{lvm}-object
#' @param links the name of the penalized links
#' @param group define the groups of the penalty
#' 
#' @examples 
#' ## no category
#' m <- lvm(Y~X1+X2+X3)
#' initGroup.lvm(m, links = c("Y~X1","Y~X2"))
#' 
#' ## categories
#' categorical(m, labels = c("A","B","C")) <- "X1"
#' initGroup.lvm(m, links = c("Y~X1","Y~X2"))
#' categorical(m, labels = c("A","B","C")) <- "X2"
#' initGroup.lvm(m, links = c("Y~X1","Y~X2"))
#' 
#' regression(m) <- Z~X1+X2+X3+eta
#' latent(m) <- ~eta
#' initGroup.lvm(m, links = c("Y~X1","Y~X2","Z~X1"))
#' 
#' 
initGroup.lvm <- function(x, links, group){
  Mlink <- getIvar.lvm(x, link = links)
  
  ## classify links according to whether or not they should be group penalized
  if(!missing(group) && identical(group, FALSE)){
    group <- NULL
    link.ordinal <- NULL
    # no group penalty
  } else if(!missing(group)){
    if(identical(group, TRUE)){
      group <- rep(1, length(links))
    }
    link.ordinal <- links[is.na(group)]
  } else { # take care of categorical variables
    link.ordinal <- links[Mlink["exogenous", ] %in% names(x$attributes$ordinalparname)]  
    group <- 1:length(link.ordinal)
  }
  link.Nordinal <- setdiff(links, link.ordinal)
  
  Mlink <- rbind(Mlink, group=NA)
  
  
  ## differentiate links and penalty according to the variable type
  if(length(link.ordinal)>0){
    Mlink["group",link.ordinal] <- group
    
    xCAT <- categorical2dummy(x, sim(x, 1))$x
    MCATlink <- getIvar.lvm(xCAT, link = setdiff(coef(xCAT), coef(x)))
    MCATlink <- rbind(MCATlink,
                      exogenous.NCAT = sapply(MCATlink["exogenous",], renameFactor, ls.levels = x$attributes$labels),
                      group = NA
    )
    
    for(iterCAT in 1:length(link.ordinal)){
      iIndex <- intersect(which(MCATlink["endogenous",] %in% Mlink["endogenous",link.ordinal[iterCAT]]),
                          which(MCATlink["exogenous.NCAT",] %in% Mlink["exogenous",link.ordinal[iterCAT]])
      )
      MCATlink["group",iIndex] <- Mlink["group",link.ordinal[iterCAT]]
    }
    
    # remove non penalized categorical variables
    MCATlink <- MCATlink[,!is.na(MCATlink["group",]),drop = FALSE]
    Mlink <- Mlink[, link.Nordinal, drop = FALSE]
  
  }else{
    MCATlink <- NULL
  }
  
  ## export
  return(list(M = Mlink,
              Mcat = MCATlink))
}

initVcoef.lvm <- function(x, link, group){
  
  
  #### create matrix
  x <- categorical2dummy(x, sim(x, 1))$x
  allCoef <- coef(x)
  n.allCoef <- length(allCoef)
  
  V <- Matrix::Matrix(rnorm(allCoef), sparse = TRUE, doDiag = FALSE,  # force to be non-symetric
                      nrow = n.allCoef, ncol = length(link),
                      dimnames = list(allCoef,NULL)
                      )
  
  #### fill matrix
  V[,] <- 0
  for(iterLink in 1:length(link)){
    V[link[iterLink],iterLink] <- as.numeric(group[iterLink])
  }
  
  return(V)
}

