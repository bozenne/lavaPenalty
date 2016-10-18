#### intialisation ####

#' @title Initialize lasso, ridge, group lasso penalty
initPenaltyL12 <- function(){
  
  penalty <- list(link = NULL,
                  group = NULL,
                  var.coef = NULL,
                  lambda1 = 0, 
                  lambda2 = 0,
                  adaptive = FALSE,
                  proxOperator = NULL,
                  objectivePenalty = NULL,
                  V = NULL)
  
  class(penalty) <- "penaltyL12"
  return(penalty)
  
}

#' @title Initialize nuclear norm penalty
initPenaltNuclear <- function(){
  
  penaltyNuclear <- list(FCTobjective = NULL,
                         FCTgradient = NULL,
                         FCThessian = NULL,
                         lambdaN  = 0, 
                         name.coef = NULL,
                         name.Y = NULL,
                         name.X = NULL,
                         nrow = NULL,
                         ncol = NULL)
  
  class(penaltyNuclear) <- "penaltyNuclear"
  return(penaltyNuclear)
  
}

#' @title Add penalty to a latent variable model
lvm2plvm <- function(x){
  
  x$penalty <- initPenaltyL12()
  x$penaltyNuclear <- initPenaltNuclear()
  class(x) <- append("plvm", class(x))
  return(x)
  
}

##' @title Penalty term for LVM
##' @aliases penalize penalize<- penalize.lvm penalize.plvm
##' @param x a lvm model
##' @param intercept should the intercept be penalized
##' @param regression should the mean parameters be penalized
##' @param variance should the variance parameters be penalized (not possible now)
##' @param lambda1 the lasso penalty
##' @param lambda2 the ridge penalty
##' @param fn_penalty user defined penalty function. Arguments coef, lambda1, lambda2.
##' @param gn_penalty first derivative of the user defined penalty function. Arguments coef, lambda1, lambda2.
##' @param hn_penalty, second derivative of the user defined penalty user defined penalty function. Arguments coef, lambda1, lambda2.
##' @param value
##' 
##' @return none
##' 
##' @examples
##' set.seed(10)
##' n <- 500
##' formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
##' lvm.modelSim <- lvm()
##' regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
##' distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
##' df.data <- sim(lvm.modelSim,n)
##' 
##' lvm.model <- lvm(formula.lvm)
##' plvm.model <- penalize(lvm.model)
##' 
##' #### regularization
##' elvm.L2 <- estimate(plvm.model,  data = df.data, lambda2 = 50)
##' elvm.L1 <- estimate(plvm.model,  data = df.data, lambda1 = 65)
##' elvm.L12 <- estimate(plvm.model,  data = df.data, lambda1 = 50, lambda2 = 50)
##' 
##' #### path regularization
##' elvm.FixedPath <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fix.sigma = TRUE,
##'                            control = list(data = df.data))
##'                            
##' elvm.L1Fixedpath <- estimate(plvm.model,  data = df.data, lambda1 = elvm.FixedPath$opt$message[3,"lambda"],
##'                         fix.sigma = TRUE)                           
##' coef(L1Fixedpath) - elvm.FixedPath$opt$message[3,-1]
##'
##' elvm.FreePath <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, control = list(data = df.data))                                                        
##' elvm.L1FreePath <- estimate(plvm.model,  data = df.data, lambda1 = elvm.FreePath$opt$message[3,"lambda"])                           
##' coef(elvm.L1FreePath) - elvm.FreePath$opt$message[3,-1]
##' 
##' @export
##' 
##' 






#### 2- Penalty functions ####

#### Lasso, group lasso, nuclear norm ####
##' @export
`penalize` <-
  function(x,...) UseMethod("penalize")

##' @export
"penalize<-" <- function (x, ..., value) {
  UseMethod("penalize<-", x)
}

##' @export
`penalize.lvm` <- function(x, value = NULL, ...){
  
  penalize(x, ...) <- value
  
  return(x)
}

`penalize<-.lvm` <- function(x,group, V, add = TRUE, reduce = FALSE,
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
    if(latent == FALSE){
      request <- paste( paste0("^",names(x$latent),"~|~",names(x$latent),"$"), collapse = "|")
      sapply("^u~", grep, x = coef(x), value = TRUE)
      ls.penaltyCoefLatent <- sapply(request, grep, x = coef(x), value = FALSE)
      index.penaltyCoef <- setdiff(index.penaltyCoef, unique(unlist(ls.penaltyCoefLatent)))
    }
    
    value <- coef(x)[index.penaltyCoef]
    
  } 
  
  #### update the name of the penalized parameters
  group.oldPenalty <- penalty(x, type = "group")
  n.oldPenalty <- length(group.oldPenalty)
  
  if(add && n.oldPenalty>0){
    
    penalty(x, type = "link") <- c(penalty(x, type = "link"), value)
    
    allGroups <- seq(0.1, 0.9, length.out = n.oldPenalty+length(value))
    allGroups[which(group.oldPenalty>=1)] <- group.oldPenalty[which(group.oldPenalty>=1)]
    penalty(x, type = "group") <- allGroups[1:n.oldPenalty]
    newGroup <-  allGroups[-(1:n.oldPenalty)]
  }else{
    penalty(x, type = "link") <- value
    
    allGroups <- seq(0.1, 0.9, length.out = length(value)) 
    penalty(x, type = "group") <- NULL
    newGroup <- allGroups
  }
  
  x$penalty$var.coef <- paste(c(endogenous(x),latent(x)),c(endogenous(x),latent(x)),sep =",") # useless
  
  #### group penalty
  if(missing(group) == FALSE){
    
    if(length(group) == length(penalty(x, type = "link"))){
      penalty(x, type = "group") <- group
    }else{
      
      if(length(group)!=1){
        stop("penalize<-.plvm: can only create one group of penalized variable at a time \n",
             "length(group): ",length(group),"\n")
      }
      
      group.id <- if(all(penalty(x, type = "group")<1)){1}else{max(penalty(x, type = "group")+1)}
      
      if(group %in% lava::vars(x)){
        newGroup[grep(group, penalty(x, type = "link"))] <- group.id
      }else{
        newGroup[] <- group.id
      }
      penalty(x, type = "group") <- c(penalty(x, type = "group"),newGroup)
    }
  }else{
    penalty(x, type = "group") <- c(penalty(x, type = "group"),newGroup)
  }
  
  #### V matrix
  if(!missing(V)){
    penalty(x, type = "V") <- V
  }else{ # erase existing settings
    V <- matrix(0, nrow = length(coef(x)), ncol = length(coef(x)))
    colnames(V) <- coef(x)
    rownames(V) <- coef(x)
    diag(V)[coef(x) %in% penalty(x, type = "link")] <- 1
    penalty(x, type = "V") <- V
  }
  
  #### penalization parameters
  if(!missing(lambda1)){
    penalty(x, type = "lambda1") <- as.numeric(lambda1)
  }
  
  if(!missing(lambda1)){
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


#' @title Remove penalty from a penalized latent variable model
#' @name cancelPenalty
#' @description Remove one or several penalties from a penalized latent variable model
#' 
#' @param x \code{plvm}-object
#' @param link the penalty that should be removed
#' @param value the penalty that should be removed
#' @param simplify if the object contain no more penalty should it be converted to a non penalized lvm
#' 
#' @examples 
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:10),y="y")
#' pm <- penalize(m)
#' 
#' cancelPenalty(pm, link = "y~x5")
#' cancelPenalty(pm) <- "y~x1"
#' cancelPenalty(pm) <- y~x1

##' @export
`cancelPenalty` <-
  function(x,...) UseMethod("cancelPenalty")

##' @export
`cancelPenalty<-` <- function (x, ..., value) {
  UseMethod("cancelPenalty<-", x)
}

`cancelPenalty.plvm` <- function(x, simplify = TRUE, link){
  cancelPenalty(x, simplify) <- link
  return(x)
}

`cancelPenalty<-.plvm` <- function(x, simplify = TRUE, value){
  
  penalty <- penalty(x, type = NULL)
  cancelPenalty(penalty) <- value
  x$penalty <- penalty
  
  if(simplify && length(penalty(x, type = "link")) == 0){
    class(x) <- setdiff(class(x), "plvm")
  }
  
  return(x)
  
}

`cancelPenalty<-.penaltyL12` <- function(x, value){
  
  link <- x$link
  index.rm <- which(link %in% value)
  
  x$group <- x$group[-index.rm]
  
  x$link <- setdiff(link, value)
  x$V[value,] <- 0
  x$V[,value] <- 0
  
  return(x)
}

#### Nuclear ####

##' @export
`penalizeNuclear` <-
  function(x,...) UseMethod("penalizeNuclear")

##' @export
"penalizeNuclear<-" <- function (x, ..., value) {
  UseMethod("penalizeNuclear<-", x)
}

##' @export
`penalizeNuclear.lvm` <- function(x, value = NULL, ...){
  
  penalizeNuclear(x, ...) <- value
  
  return(x)
}

##' @export
`penalizeNuclear.plvm` <- `penalizeNuclear.lvm`

##' @export
`penalizeNuclear<-.lvm` <- function(x, ..., value){
  
  ## convert to plvm
  x <- lvm2plvm(x)
  
  ## main
  penalizeNuclear(x, ...) <- value
  
  ## export
  return(x)
}

##' @export
`penalizeNuclear<-.plvm` <- function(x,  coords, objective = NULL, gradient = NULL, hessian = NULL, lambda = NULL, ..., value){
  
  ## add objective/gradient
  if(is.null(objective) && is.null(x$penaltyNuclear$FCTobjective)){
    x$penaltyNuclear$FCTobjective <- lvGaussian
  }
  
  if(is.null(gradient) && is.null(x$penaltyNuclear$FCTgradient)){
    x$penaltyNuclear$FCTgradient <- scoreGaussian
  }
  
  if(is.null(hessian) && is.null(x$penaltyNuclear$FCThessian)){
    x$penaltyNuclear$FCThessian <- hessianGaussian
  }
  
  if(!is.null(objective)){
    x$penaltyNuclear$FCTobjective <- objective
  }
  if(!is.null(gradient)){
    x$penaltyNuclear$FCTgradient <- gradient
  }
  if(!is.null(hessian)){
    x$penaltyNuclear$FCThessian <- hessian
  }
  
  ## lambda
  if(!is.null(lambda)){
    x$penaltyNuclear$lambdaN <- lambda  
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
  x$penaltyNuclear$name.coef <- paste0(name.Y,"~",name.X)
  x$penaltyNuclear$name.Y <- name.Y
  x$penaltyNuclear$name.X <- name.X
  x$penaltyNuclear$ncol <- ncol  
  x$penaltyNuclear$nrow <- nrow  
  
  #### export
  return(x)
}

