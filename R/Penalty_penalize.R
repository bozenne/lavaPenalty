#### intialisation ####

# {{{ initPenalty

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
                    lambdaG = 0,
                    adaptive = FALSE,
                    VelasticNet = NULL,
                    Vgroup = NULL)
  
  class(penalty) <- "penaltyL12"
  return(penalty)
  
}

#' @rdname initPenalty
initPenaltNuclear <- function(){
  
    penaltyNuclear <- list(lambdaN  = 0,
                           nrow = NULL,
                           ncol = NULL,                           
                           name.reduce = NULL,
                           endogeneous = NULL,
                           link = NULL)
  
  class(penaltyNuclear) <- "penaltyNuclear"
  return(penaltyNuclear)
  
}

# }}}

#### convertion ####

# {{{ lvm2plvm

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

# }}}

#### specification ####

# {{{ penalize

# {{{ doc
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
#' @param lambdaG group lasso penalization parameter
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
#' @export
`penalize` <-
  function(x,...) UseMethod("penalize")
# }}}


#' @rdname penalize
#' @export
"penalize<-" <- function (x, ..., value) {
  UseMethod("penalize<-", x)
}

#' @rdname penalize
#' @export
`penalize.lvm` <- function(x, value = NULL, ...){
  
  penalize(x, ...) <- value
  
  return(x)
}

#' @rdname penalize
#' @export
`penalize<-.lvm` <- function(x, group, V, add = TRUE, reduce = FALSE,
                             lambda1, lambda2, lambdaG, lambdaN, adaptive,
                             intercept = FALSE, regression = TRUE, variance = FALSE, covariance = FALSE, latent = FALSE,
                             value){
  
  ## convert to plvm
  if("plvm" %in% class(x) == FALSE){
    x <- lvm2plvm(x)
  }
  
    #### find coefficients from value
  if(!is.null(value)){
    
    if("formula" %in% class(value)){
        value <- unlist(lapply(lava.reduce::initVar_link(value, format = "formula"), formula2character))
    }else if(is.list(value)){
      value <- unlist(lapply(value, function(v){
        if("formula" %in% class(v)){unlist(lapply(lava.reduce::initVar_link(v, format = "formula"), formula2character))}else{v}
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
      index.penaltyCoef <- c(index.penaltyCoef, coefIntercept(x))  
    }
    if(regression == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, coefReg(x)) 
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
  } else {
     stop("value must contains the names of the penalized associations \n")
  }
  
    #### update the V matrices
    if(!missing(V)){ ## predefined V
        if("Matrix" %in% is(V)){
            stop("V must herit from the class Matrix \n")
        }
        penalty(x, type = "V") <- V
    }else{ ## group penalty
        resInit <- initGroup.lvm(x, links = value, group = group)

        if(!is.null(resInit$Mcat)){
            Vgroup.old <- penalty(x, type = "Vgroup")
            if(add && !is.null(Vgroup.old)){
                newV <- initVcoef.lvm(x,
                                      link = rownames(resInit$Mcat),
                                      group = resInit$Mcat[,"group"])
                newV@x <- newV@x + max(Vgroup.old@x)
                penalty(x, type = "Vgroup") <- cbind(Vgroup.old,newV)
        
      }else{
        penalty(x, type = "Vgroup") <- initVcoef.lvm(x,
                                                     link = rownames(resInit$Mcat),
                                                     group = resInit$Mcat[,"group"])
      }
    }

        if(NCOL(resInit$M)>0){  ## elastic net
            penalty(x, type = "VelasticNet") <- initVcoef.lvm(x,
                                                              link = rownames(resInit$M),
                                                              group = 1:NROW(resInit$M))
        }
    
  }
  
    #### penalization parameters
    if(!missing(lambda1)){
        penalty(x, type = "lambda1") <- as.numeric(lambda1)
    }
  
    if(!missing(lambda2)){
        penalty(x, type = "lambda2") <- as.numeric(lambda2)
    }

    if(!missing(lambdaG)){
        penalty(x, type = "lambdaG") <- as.numeric(lambdaG)
    }

    if(!missing(lambdaN)){
        penalty(x, type = "lambdaN") <- as.numeric(lambdaN)
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

# }}}

# {{{ penalizeNuclear

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
`penalizeNuclear<-.lvm` <- function(x, coords = NULL, lambdaN = NULL, ..., value){
    symbols <- lava.options()$Nuclear$symbols
    
    #### read and reshape arguments
    # convert to plvm
    if("plvm" %in% class(x) == FALSE){
        x <- lvm2plvm(x)
    }

    # lambda
    if(!is.null(lambdaN)){
        penalty(x, type = "lambdaN", nuclear = TRUE) <- lambdaN
    }

    #### value argument
    # can be a matrix: extract coefficients names and create coords
    if(is.null(coords)){
        if(!is.matrix(value)){
            stop("When argument \'coords\' is missing argument \'value\' must be a matrix \n")
        }
        coords <- which(value != 0, arr.ind = TRUE)
        value <- as.vector(value)
    }

    # can be a vector containing the coefficient names: convert to a unique formula
    if(length(value) > 1 && ("character" %in% class(value))){
        value <- lava.reduce::combine.formula(value)
        if(length(value)>1){
            stop("value must correspond to only one outcome \n")
        }
        value <- value[[1]]
    }

    # can be a formula
    if("formula" %in% class(value) == FALSE){
        stop("If not a matrix, value must be a formula or a vector containing the name of the coefficients \n")
    }
    resTempo <- lava.reduce::initVar_link(value, repVar1 = TRUE)
    name.Y <- resTempo$var1
    name.X <- resTempo$var2

    # consistency with coords
    if(NROW(coords) != length(name.X)){
        stop("Inconsistency between argument \'coords\' and argument \'value\' \n",
             "NROW(coords): ", NROW(coords),"\n",
             "length(value): ",length(name.X),"\n")
    }
    ncol <- max(coords[,1])
    nrow <- max(coords[,2])  
    if(ncol*nrow != NROW(coords)){
        stop("argument \'coords\' does not corresponds to a full matrix \n")
    }
    if(any(duplicated(coords)) ){
        stop("argument \'coords\' must not contain duplicated coordinates \n")
    }

    # consistency with the LVM
    if(unique(name.Y) %in% endogenous(x) == FALSE){
        stop("penaltyNuclear: the dependent variable in formula must already be in the model \n")
    }

    if(any(name.X %in% vars(x) == TRUE)){
        stop("penaltyNuclear: the independent variable in formula must not be in the model \n",
             "existing variables: ",paste(name.X[name.X %in% vars(x) == TRUE],collapse = " "),"\n")
    }

    #### Define penalty
    name.reduce <- paste0(symbols[1],LCSseq(name.X),symbols[2])
    vec.coef <- rep("NA", NROW(coords))
    vec.coef[coords[,1]+ncol*(coords[,2]-1)] <- paste0(name.Y,lava.options()$symbols[1], name.X)
    
    #### update penalty    
    penalty(x, type = "nrow", nuclear = TRUE) <- c(penalty(x, type = "nrow", nuclear = TRUE),
                                                   nrow)
    penalty(x, type = "ncol", nuclear = TRUE) <- c(penalty(x, type = "ncol", nuclear = TRUE),
                                                   ncol)
    penalty(x, type = "link", nuclear = TRUE) <- c(penalty(x, type = "link", nuclear = TRUE),
                                                    list(vec.coef))
    penalty(x, type = "name.reduce", nuclear = TRUE) <- c(penalty(x, type = "name.reduce", nuclear = TRUE),
                                                          name.reduce)
    penalty(x, type = "endogeneous", nuclear = TRUE) <- c(penalty(x, type = "endogeneous", nuclear = TRUE),
                                                          name.Y[1])
  
    #### update object
    x <- lava.reduce:::lvm2reduce(x)
    x <- regression(x, from = name.X, to = name.Y[1],
                    reduce = name.reduce) 
  
  #### export
  return(x)
}

# }}}


#### init ####

# {{{ initGroup.lvm

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
initGroup.lvm <- function(x, links, data = sim(x, 1), group){
    Mlink <- getIvar.lvm(x, link = links, format = "data.frame")
    ## classify links according to whether or not they should be group penalized
    if(missing(group)){
        link.ordinal <- links[Mlink[,"exogenous"] %in% names(x$attributes$ordinalparname)]
        group <- 1:length(link.ordinal)
    }else if(identical(group, FALSE) ){
        group <- NULL
        link.ordinal <- NULL
        # no group penalty
    } else if(!missing(group)){
        if(identical(group, TRUE)){
            group <- rep(1, length(links))
        }
        link.ordinal <- links[is.na(group)]
    }
    link.Nordinal <- setdiff(links, link.ordinal)
  
    Mlink <- cbind(Mlink, group=NA)

    
    
  ## differentiate links and penalty according to the variable type
  if(length(link.ordinal)>0){
     Mlink[link.ordinal,"group"] <- group
     xCAT <- lava_categorical2dummy(x, data)$x
     MCATlink <- getIvar.lvm(xCAT, format = "data.frame")
     MCATlink[setdiff(coef(xCAT),coef(x)),"type"] <- "categorical"
     MCATlink <- MCATlink[MCATlink[,"type"] == "categorical",]
     MCATlink <- cbind(MCATlink,
                       group = as.numeric(NA),
                       originalLink = as.character(NA),
                       stringsAsFactors = FALSE)
     
    for(iterCAT in 1:length(link.ordinal)){
      
        endoCAT <- Mlink[link.ordinal[iterCAT],"endogenous"]
        exoCAT <- Mlink[link.ordinal[iterCAT],"exogenous"]
        if(exoCAT %in% names(x$attributes$labels)){
            exoCAT <- paste0(exoCAT, x$attributes$labels[[exoCAT]])
        }
      
        iIndex <- intersect(which(MCATlink[,"endogenous"] %in% endoCAT),
                            which(MCATlink[,"exogenous"] %in% exoCAT)
                            )
        MCATlink[iIndex,"group"] <- Mlink[link.ordinal[iterCAT],"group"]
        MCATlink[iIndex,"originalLink"] <- link.ordinal[iterCAT]
        
    }
     
    # remove non penalized categorical variables
    MCATlink <- MCATlink[!is.na(MCATlink[,"group"]),,drop = FALSE]
    Mlink <- Mlink[link.Nordinal,, drop = FALSE]
  
  }else{
    MCATlink <- NULL
  }

  ## export
  return(list(M = Mlink,
              Mcat = MCATlink))
}

# }}}

# {{{ initVcoef.lvm
initVcoef.lvm <- function(x, link, group){
  
  
    #### create matrix
    x <- lava_categorical2dummy(x, sim(x, 1))$x
    allCoef <- coef(x)
    n.allCoef <- length(allCoef)
    n.groups <- length(unique(group))

    V <- Matrix::Matrix(rnorm(allCoef), sparse = TRUE, doDiag = FALSE,  # force to be non-symetric
                        nrow = n.allCoef, ncol = n.groups,
                        dimnames = list(allCoef,NULL)
                        )

    #### fill matrix
    V[,] <- 0
    for(iterLink in 1:length(link)){ # iterLink <- 1
        V[link[iterLink],group[iterLink]] <- 1#as.numeric(group[iterLink])
    }
 
  return(V)
}
# }}}

