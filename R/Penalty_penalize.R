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
  
  penaltyNuclear <- list(link = NULL,
                         lambdaN  = 0, 
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
#' #### regularization
#' elvm.L2 <- estimate(plvm.model,  data = df.data, lambda2 = 50)
#' elvm.L1 <- estimate(plvm.model,  data = df.data, lambda1 = 65)
#' elvm.L12 <- estimate(plvm.model,  data = df.data, lambda1 = 50, lambda2 = 50)
#' 
#' #### path regularization
#' elvm.FixedPath <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fix.sigma = TRUE,
#'                            control = list(data = df.data))
#'                            
#' elvm.L1Fixedpath <- estimate(plvm.model,  data = df.data, lambda1 = elvm.FixedPath$opt$message[3,"lambda"],
#'                         fix.sigma = TRUE)                           
#' coef(L1Fixedpath) - elvm.FixedPath$opt$message[3,-1]
#'
#' elvm.FreePath <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, control = list(data = df.data))                                                        
#' elvm.L1FreePath <- estimate(plvm.model,  data = df.data, lambda1 = elvm.FreePath$opt$message[3,"lambda"])                           
#' coef(elvm.L1FreePath) - elvm.FreePath$opt$message[3,-1]
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
      request <- paste( paste0("^",names(x$latent),"~|~",names(x$latent),"$"), collapse = "|")
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

#' @export
`penalizeNuclear.plvm` <- `penalizeNuclear.lvm`

#' @export
`penalizeNuclear<-.lvm` <- function(x, ..., value){
  
  ## convert to plvm
  x <- lvm2plvm(x)
  
  ## main
  penalizeNuclear(x, ...) <- value
  
  ## export
  return(x)
}

#' @name penalize
#' @export
`penalizeNuclear<-.plvm` <- function(x, coords, lambdaN = NULL, ..., value){
  
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
#' @export
`cancelPenalty` <-
  function(x,...) UseMethod("cancelPenalty")

#' @export
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

# `cancelPenalty<-.penaltyNuclear` TODO!!!!

