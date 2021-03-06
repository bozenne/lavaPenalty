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
  
    penalty <- list(lambda1 = numeric(0),
                    lambda2 = numeric(0),
                    lambdaG = numeric(0),
                    Vlasso = NULL,
                    Vridge = NULL,
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
#' @description Conversion to a penalized latent variable model
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
#' @param coords the (spatial) position of the links to penalize
#' @param Vlasso the matrix defining lasso penalties. Otherwise must be logical to indicate whether to add lasso penalties.
#' @param Vridge the matrix defining ridge penalties. Otherwise must be logical to indicate whether to add ridge penalties.
#' @param Vgroup the matrix defining group lasso penalties. Otherwise must be logical to indicate whether to add group lasso penalties.
#' @param add should value be added to the existing penalty term ? Otherwise it will overwrite it.
#' @param reduce should for each regression the penalised link be aggregated into a linear predictor.
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

#' @rdname penalize
#' @export
"penalize<-" <- function (x, ..., value) {
  UseMethod("penalize<-", x)
}
# }}}

# {{{ penalize.lvm

#' @rdname penalize
#' @export
`penalize.lvm` <- function(x, value = NULL, ...){
  
  penalize(x, ...) <- value
  
  return(x)
}

# }}}
# {{{ penalize<-.lvm
#' @rdname penalize
#' @export
`penalize<-.lvm` <- function(x, group, add = TRUE, reduce = FALSE,
                             Vlasso = TRUE, Vridge = TRUE, Vgroup = TRUE, 
                             intercept = FALSE, regression = TRUE, variance = FALSE, covariance = FALSE, latent = FALSE,
                             value, ...){

    ## convert to plvm
    if("plvm" %in% class(x) == FALSE){
        x <- lvm2plvm(x)
    }

    test.V <- !is.logical(Vlasso)  || !is.logical(Vridge) || !is.logical(Vgroup)
    if(!is.null(value) && test.V){
        stop("argument \'value\' must be NULL when arguments \'Vlasso\', \'Vridge\', or \'Vgroup\' are matrices \n")
    }

    #### update of the matrix
    if(test.V){
        if(!missing(Vlasso)){ ## predefined V
            if("Matrix" %in% is(Vlasso)){
                stop("Vlasso must herit from the class Matrix \n")
            }
            penalty(x, type = "Vlasso") <- Vlasso
        }
        if(!missing(Vridge)){ ## predefined V
            if("Matrix" %in% is(Vridge)){
                stop("Vridge must herit from the class Matrix \n")
            }
            penalty(x, type = "Vridge") <- Vridge
        }
        if(!missing(Vgroup)){ ## predefined V
            if("Matrix" %in% is(Vgroup)){
                stop("Vgroup must herit from the class Matrix \n")
            }
            penalty(x, type = "Vgroup") <- Vgroup
        }
    } else {
        #### otherwise find coefficients from value
        if(!is.null(value)){

            value <- lavaReduce::initVar_links(value, format = "txt.formula")    
    
            if(any(value %in% coef(x) == FALSE)){
                stop("penalty<-.lvm: coefficients to be penalized do not match those of the model\n",
                     "unknown coefficients: ",paste(value[value %in% coef(x) == FALSE], collapse = " "),"\n",
                     "available coefficients: ",paste(coef(x)[coef(x) %in% value == FALSE], collapse = " "),"\n")
            }
    
        } else {
    
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
        }

        ## group penalty
        resInit <- initGroup.lvm(x, link = value, group = group)
        if(Vgroup && resInit[type == "categorical",.N]>0){
            newV <- initVcoef.lvm(x,
                                  link = resInit[type == "categorical"][["link"]],
                                  group = resInit[type == "categorical"][["group"]])
            penalty(x, type = "Vgroup", add = add) <- newV        
        }

        if(resInit[type == "continuous",.N]>0){  ## elastic net
            # lasso
            if(Vlasso){
            Vlasso <- initVcoef.lvm(x,
                                    link = resInit[type == "continuous"][["link"]],
                                    group = 1:resInit[type == "continuous",.N])           
            penalty(x, type = "Vlasso", add = add) <- Vlasso
            }
            
            # ridge
            if(Vridge){
            Vridge <- initVcoef.lvm(x,
                                    link = resInit[["link"]],
                                    group = 1:resInit[,.N])            
            penalty(x, type = "Vridge", add = add) <- Vridge
            }
        }
    
    }
    
    #### reduce 
    if(reduce){
        x <- reduce(x)
    }
  
  
    #### export
    return(x)
}
# }}}

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
        value <- lavaReduce::combine.formula(value)
        if(length(value)>1){
            stop("value must correspond to only one outcome \n")
        }
        value <- value[[1]]
    }

    # can be a formula
    if("formula" %in% class(value) == FALSE){
        stop("If not a matrix, value must be a formula or a vector containing the name of the coefficients \n")
    }
    resTempo <- lavaReduce::initVar_link(value, repVar1 = TRUE)
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
    x <- lavaReduce::lvm2reduce(x)
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
#' dt <- initGroup.lvm(m, links = c("Y~X1","Y~X2"))
#' dt
#' 
#' ## categories
#' categorical(m, labels = c("A","B","C")) <- "X1"
#' dt <- initGroup.lvm(m, links = c("Y~X1","Y~X2"))
#' dt
#' categorical(m, K=2, labels = 1:2) <- "X2"
#' dt <- initGroup.lvm(m, links = c("Y~X1","Y~X2"))
#' dt
#' 
#' regression(m) <- Z~X1+X3+X5+eta
#' latent(m) <- ~eta
#' dt <- initGroup.lvm(m, links = c("Y~X1","Y~X2","Z~X1"))
#' dt
initGroup.lvm <- function(x, links, group){

    dt.link <- getIvar.lvm(x, link = links, format = "data.table")
    dt.link[, group := as.numeric(NA)]

    ## update link according to categorical variables
    newlink <- dt.link[["link"]]

    ## classify links according to whether or not they should be group penalized
    if(missing(group)){
        dt.link[type == "categorical", group := as.numeric(.GRP), by = c("endogenous","exogenous")]
    }else if(identical(group, FALSE) ){
        # no group penalty
    }else if(identical(group, TRUE)){
        dt.link[type == "categorical", group := as.numeric(1)]
    }else {
        if(length(group)==1){group <- rep(group, length(links))}
        groupClean <- group[!is.na(group)] # necessary because otherwise ambiguity dt var
        linksClean <- links[!is.na(group)]
        dt.link[match(linksClean,link), group := groupClean]        
    }    

    ## export
    return(dt.link)
}

# }}}

# {{{ initVcoef.lvm
initVcoef.lvm <- function(x, link, group){
    
    #### create matrix
    x <- lava_categorical2dummy(x, sim(x, 1))$x
    allCoef <- coef(x)
    n.allCoef <- length(allCoef)
    n.groups <- length(unique(group))

    V <- Matrix::Matrix(rnorm(n.allCoef), sparse = TRUE, doDiag = FALSE,  # force to be non-symetric
                        nrow = n.allCoef, ncol = n.groups,
                        dimnames = list(allCoef,NULL)
                        )

    #### fill matrix
    V[] <- 0
    for(iterLink in 1:length(link)){ # iterLink <- 1
        V[link[iterLink],group[iterLink]] <- 1#as.numeric(group[iterLink])
    }
 
  return(V)
}
# }}}

# {{{ initVcoef.multigroup
initVcoef.multigroup <- function(x, type){
    n.models <- length(x$lvm)

    ##  coef names new
    allCoef <- x$name
    n.allCoef <- length(allCoef)

    x.coef <- mapply(c, x$meanlist, x$parlist, SIMPLIFY = FALSE)
    posTempo <- unlist(lapply(strsplit(allCoef,split="@"),"[",1))
    for(iModel in 1:n.models){ # iModel <- 1
        x.coef[[iModel]] <- setNames(allCoef[posTempo==as.character(iModel)],
                                     names(x.coef[[iModel]]))
    }

    ## coef names old
    ls.coef <- lapply(x$lvm,coef) # list of lvm
    ls.penalty <- lapply(x$lvm, penalty)
    max.penalty <- lapply(ls.penalty,function(df){max(df$group,na.rm=TRUE)})
    n.penalty <- cumsum(c(0,unlist(max.penalty)))

    ## find penalty
    ls.newPenalty <- lapply(1:n.models, function(iModel){ # iModel <- 2
        indexCoef.tempo <- match(ls.penalty[[iModel]][penalty==type,link],ls.coef[[iModel]])
        name.tempo <- names(ls.coef[[iModel]])[indexCoef.tempo]

        list(name = x.coef[[iModel]][name.tempo],
             group = ls.penalty[[iModel]][penalty==type,group] + n.penalty[iModel])
    })
    
    ## remove duplicated penalties
    allPenalty <- unlist(lapply(ls.newPenalty,"[[","name"))
    allGroup <- unlist(lapply(ls.newPenalty,"[[","group"))
    index.keep <- which(!duplicated(allPenalty))
    newPenalty <- allPenalty[index.keep]
    group <- allGroup[index.keep]
  
    ## create V
    V <- Matrix::Matrix(rnorm(n.allCoef), sparse = TRUE, doDiag = FALSE,  # force to be non-symetric
                        nrow = n.allCoef, ncol = max(group),
                        dimnames = list(allCoef,NULL)
                        )
  
    #### fill matrix
    V[] <- 0
    for(iterLink in 1:length(newPenalty)){ # iterLink <- 1
        V[newPenalty[iterLink],group[iterLink]] <- 1#as.numeric(group[iterLink])
    }
    return(V)
}
# }}}
