### Penalty_initializer.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar  6 2017 (11:51) 
## Version: 
## last-updated: apr 20 2017 (15:34) 
##           By: Brice Ozenne
##     Update #: 242
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ initialize

# {{{ doc
#' @title Intialization functions for plvm models
#' @name initialize
#' @description
#' initializeFactor.penaltyL12: update the lasso penalty to group lasso penalty in presence of categorical variables.
#' initialize.penaltyL12: estimate starting value for the plvm.
#' initialize.penaltyNuclear: estimate starting value for the plvm.
#' initialize.start: estimate starting value for the plvm.
#' 
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param name.coef the name of the parameters
#' @param trace should the user be told that some penalties have been updated ?
#'
#' @details
#' Remove penalties corresponding to reference links
#' 
#' @examples
#'
#' #### initialize.data ####
#'
#'
#' #### initializeFactor.penaltyL12 ###
#' m <- lvm(Y~X1+X2+X3+X4)
#' pm <- penalize(m)
#' initializeFactor.penaltyL12(pm, data = sim(m, 1), trace = TRUE)
#' 
#' mCAT <- m
#' categorical(mCAT, K = 3, labels = letters[1:3]) <- "X1"
#' categorical(mCAT, K = 2, labels = letters[1:2]) <- "X2"
#'
#' pm <- penalize(m)
#' initializeFactor.penaltyL12(pm, data = sim(mCAT, 1), trace = TRUE)
#'
#' pmCAT <- penalize(mCAT)
#' initializeFactor.penaltyL12(pmCAT, data = sim(mCAT, 1), trace = TRUE)
#' 
#' #### initialize.penaltyL12 #### 
#' m <- lvm()
#' regression(m) <- y1~x1+x2+x3+x4
#'
#' #### elastic net
#' e <- estimate(m, sim(m, 1e2))
#' pm <- penalize(m, lambda1 = 2, lambda2 = 1.5)
#' pen12 <- initialize.penaltyL12(pm$penalty, name.coef = names(coef(e)))
#' pen12
#'
#' i <- initializer.plvm(pm, data = sim(m, 1e2),
#'                  regularizationPath = FALSE, constrain.variance = FALSE)
#' 
#' #### group lasso
#' categorical(m, labels = c("A","B","C")) <- "x1"
#' categorical(m, labels = c("A","B","C")) <- "x2"
#'
#' e <- estimate(m, sim(m, 1e2))
#' pm <- penalize(m, lambdaG = 5)
#' pen12 <- initializer.penaltyL12(pm$penalty, name.coef = names(coef(e)))
#' pen12
#' estimate(pm, data = sim(pm, 1e2))
#'
#' i <- initializer.plvm(pm, data = sim(m, 1e2),
#'                  regularizationPath = FALSE, constrain.variance = FALSE)
#' 
#' # nuclear norm
#' m <- lvm()
#' coords <- expand.grid(x = 1:10, y = 1:10)
#' m <- regression(m, y = "y1", x = paste0("z",1:100))
#'
#' mNuclear <- lvm(y1 ~ x1 + x2)
#' penalizeNuclear.penaltyL12(mNuclear, coords = coords, lambdaN = 10) <- coefReg(m, value = TRUE)
#' mNuclear
#'
#' penN <- initializer.plvm(mNuclear$penaltyNuclear, name.coef = coef(mNuclear))
#' penN
#' 
# }}}

# {{{ initializeFactor.penaltyL12
#' @rdname initialize
initializeFactor.penaltyL12 <- function(x, data, trace){

    table.penalty <- penalty(x, nuclear = FALSE) 
    var.penalty <- unique(table.penalty$exogenous)        
    test.factor <- sapply(var.penalty, function(v){is.factor(data[[v]])})
    
    ## search for links that should be moved from lasso to group lasso
    if(any(test.factor)){
       
        var.factor <- unique(var.penalty[test.factor])
        table.penalty.factor <- table.penalty[table.penalty$exogenous %in% var.factor & penalty == "lasso"]

        ## remove lasso link
        cancelPenalty(x, rm.lasso = TRUE, rm.ridge = FALSE, rm.group = FALSE) <- table.penalty.factor[["link"]]

        ## update the lvm with the categorical variables
        for(f in var.factor){ # f <- "X1"
            if(trace>0){cat("convert the lasso penalty on the categorical variable ",f," to a group penalty \n")}
            categorical(x, labels = levels(as.factor(data[[f]]))) <- as.formula(paste0("~",f))
        }
            
        ## update the penalty according to the categorical variables
        resInit <- initGroup.lvm(x, links = table.penalty.factor$link)
        # same as `penalize<-.lvm`
        newV <- initVcoef.lvm(x,
                              link = resInit[type == "categorical"][["link"]],
                              group = resInit[type == "categorical"][["group"]])
        penalty(x, type = "Vgroup", add = TRUE) <- newV
    }

    #### export
    return(x)
}
# }}}

# {{{ initialize.penaltyL12

#' @rdname initialize
initialize.penaltyL12 <- function(x, name.x, model.x,
                                  regularizationPath){
  
    n.coef <- length(name.x)
    table.penalty <- penalty(x)
    table.penalty[, model := model.x[link]]
    
    ## check penalized links corresponding to existing parameters
    if(any(table.penalty$link %in% name.x == FALSE)){
        rm.penalty <- setdiff(table.penalty$link, name.x)
        stop("Penalized links not found in coef  \n",
             "links: ",paste(rm.penalty, collapse = " "),"\n",
             "Could be because parameters are used as a reference")
    }
    
    ## remove extra parameters in penalty
    # since some parameter may be fixed as reference
    for(iType in c("Vlasso","Vridge","Vgroup")){ # iType <- "Vlasso"
        V <- penalty(x, type = iType)[[iType]]        
        if(!is.null(V) && any(rownames(V) %in% name.x == FALSE)){
            index.rm <- which(rownames(V) %in% name.x == FALSE)
            penalty(x, type = iType, add = FALSE) <- V[-index.rm,]
        }
    }

    if(regularizationPath){ ## create matrix

        ## extract V matrix
        Vlasso <- penalty(x, type = "Vlasso")$Vlasso                
        if(table.penalty[penalty == "ridge",.N]>0){
            index.penalty2 <- match(table.penalty[penalty=="ridge",link], name.x)
            lambda2 <- table.penalty[penalty=="ridge",lambda]
            index.model2 <- table.penalty[penalty=="ridge",model]
        }else{
            index.penalty2 <- NULL
            index.model2 <- NULL
            lambda2 <- NULL
        }

        res <- list(Vlasso = Vlasso,
                    lambda2 = lambda2, index.penalty2 = index.penalty2, index.model2 = index.model2)
        
    }else{ ## locate the penalized parameters

        ## lasso
        if(table.penalty[penalty == "lasso",.N]>0){
            index.penalty1 <- match(table.penalty[penalty=="lasso",link], name.x)
            lambda1 <- table.penalty[penalty=="lasso",lambda]
            index.model1 <- table.penalty[penalty=="lasso",model]
        }else{
            index.penalty1 <- NULL
            index.model1 <- NULL
            lambda1 <- NULL
        }

        ## ridge
        if(table.penalty[penalty == "ridge",.N]>0){
            index.penalty2 <- match(table.penalty[penalty=="ridge",link], name.x)
            lambda2 <- table.penalty[penalty=="ridge",lambda]
            index.model2 <- table.penalty[penalty=="ridge",model]
        }else{
            index.penalty2 <- NULL
            index.model2 <- NULL
            lambda2 <- NULL
        }
    
        ## define groups of parameters for the group lasso
        if(table.penalty[penalty == "group lasso",.N]>0){            
            index.penaltyG  <- table.penalty[penalty == "group lasso",.(.(setNames(match(link,name.x),link))), by  = group]$V1
            # same lambda for all parameters in the same group: no need for a list
            lambdaG  <- table.penalty[penalty == "group lasso",lambda[1], by  = group]$V1
            index.modelG <- table.penalty[penalty=="group lasso",model]
        }else{
            index.penaltyG <- NULL
            index.modelG <- NULL
            lambdaG <- NULL
        }

        res <- list(lambda1 = lambda1, index.penalty1 = index.penalty1, index.model1 = index.model1,
                    lambda2 = lambda2, index.penalty2 = index.penalty2, index.model2 = index.model2,
                    lambdaG = lambdaG, index.penaltyG = index.penaltyG, index.modelG = index.modelG)
    }
    return(res)
}

# }}}

# {{{ initialize.start
#' @rdname initialize
initialize.start <- function(x, data, regularizationPath, increasing,
                             constrain.variance, name.variance,
                             ...){

    data <- data[,manifest(x, lp = FALSE, xlp = TRUE), drop = FALSE]
  
    names.coef <- coef(x)
    # names.coef  <- coef(lava::fixsome(x, measurement.fix = TRUE, n = NROW(data), debug = FALSE))
    n.coef <- length(names.coef)
    n.data <- NROW(data)

    # {{{ compute start - low dimensional (keep all links unpenalized)

    test.ridge <- length(penalty(x, "lambda2")$lambda2)>0

    if(regularizationPath && increasing){
        if(test.ridge){
            
        suppressWarnings(
            initLVM <- tryCatch(coef(estimate(x = x, data = data,
                                         regularizationPath = FALSE,
                                         lambda1 = 0,
                                         lambda2 = penalty(x, "lambda2")$lambda2,
                                         ...),
                                error = function(e){e}
                                ))
        )
        } else {
            x0 <- x
            class(x0) <- setdiff(class(x0),"plvm")
        
            suppressWarnings(
                initLVM <- tryCatch(estimate(x = x0, data = data, quick = TRUE, ...),
                                    error = function(e){e}
                                    )
            )
        }

        if(all("error" %in% class(initLVM) == FALSE)){        
            start <- initLVM                 
        }else{ 
            start <- NULL
        }
    }

    # }}}

    # {{{ compute start - high dimensional (remove penalized links, all lambda are infinity)

    if(!regularizationPath || !increasing){
        start <- setNames(rep(0,n.coef),names.coef)

        ## simplify the LVM
        x0 <- x
        class(x0) <- setdiff(class(x0),"plvm")
    
        # remove penalized regression links
        link.penalty <- unique(penalty(x, nuclear = FALSE)$link)
        linkReg.penalty <- intersect(coefReg(x, value = TRUE),
                                     link.penalty)
        if(!is.null(linkReg.penalty)){
            ls.formula <- lavaReduce::combine.formula(linkReg.penalty)
            for(iter_link in ls.formula){
                x0 <- cancel(x0, iter_link)
            }
        }
    
        # nuclearNorm
        reg.penalty <- intersect(coefReg(x, value = TRUE),
                                 penalty(x, type = "link", nuclear= TRUE))
        if(!is.null(reg.penalty)){
            ls.formula <- lavaReduce::combine.formula(reg.penalty)
            for(iter_link in ls.formula){
                x0 <- cancel(x0, iter_link)
            }
        }
    
        # remove penalised LV [TODO]
        x0 <- clean(x0, rm.endo = FALSE)
    
        ## estimate the model
        suppressWarnings(
            newCoef <- estimate(x = x0, data = data, quick = TRUE, ...)
        )
        newCoef <- newCoef[names(newCoef) %in% names(start)] # only keep relevant parameters
        start[names(newCoef)] <- newCoef
    }

    # }}}
 
    # export
    return(start)
}

# }}}

# {{{ initialize.penaltyNuclear
#' @rdname initialize
initializer.penaltyNuclear <- function(x, name.coef){
    n.coef <- length(name.coef)
    
    name.allPenalty <- penalty(x, type = "link")$link
    lambdaN <- penalty(x, type = "lambdaN")$lambdaN
    if(length(name.allPenalty)==0){
        return(list(lambdaN = NULL, ncol = NULL, nrow = NULL, index.penaltyN = NULL))
    }

    n.nuclear <- length(name.allPenalty)
    vec.lambdaN <- setNames(rep(0, n.coef), name.coef)
    vec.ncol <- penalty(x, type = "ncol")
    vec.nrow <- penalty(x, type = "nrow")
    index.penaltyN <- vector(mode = "list", length= n.nuclear)
    
    for(iNuclear in 1:n.nuclear){

        if(any(name.allPenalty[[iNuclear]] %in% name.coef == FALSE)){
            stop("initPenalty: some penalty will not be applied because the corresponding parameter is used as a reference \n",
                 "non-applied penalty: ",paste(setdiff(name.allPenalty[[iNuclear]], name.coef), collapse = " "),"\n")
        }
        vec.lambdaN[name.allPenalty[[iNuclear]]] <- lambdaN[iNuclear]
        index.penaltyN[[iNuclear]] <- setNames(match(name.allPenalty[[iNuclear]],
                                                     name.coef), name.allPenalty[[iNuclear]])

    
    }
        
    return(list(lambdaN = vec.lambdaN, ncol = vec.ncol, nrow = vec.nrow,
                index.penaltyN = index.penaltyN))
}
# }}}

# }}}


# {{{ initialize.data
#' @title Scale the dataset and convert categorical variables to indicators
#' @description Scale the dataset and convert categorical variables to indicators
#' 
#' @name initializeData
#' 
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param method.center function used to center the data
#' @param method.scale function used to scale the data

`initializeData` <-
  function(x,...) UseMethod("initializeData")


#' @rdname initializeData
initializeData.plvm <- function(x, data, method.center = "mean", method.scale = "sd"){

    if(any(manifest(x, lp = FALSE) %in% names(data) == FALSE)){
        stop("prepareData.lvm: arguments \'data\' and \'x\' are incompatible \n",
             "variables: ",paste(manifest(x)[manifest(x) %in% names(data) == FALSE], collapse = " ")," not found in \'data\' \n")
    }

    #### convert categorical variables to dummy variables
    resC2D <- lava_categorical2dummy(x, data)
  
    index.numeric <- intersect(manifest(x, lp = FALSE), manifest(resC2D$x, lp = FALSE))
    indexOld.factor <- setdiff(manifest(x, lp = FALSE),  manifest(resC2D$x, lp = FALSE))
    indexNew.factor <- setdiff(manifest(resC2D$x, lp = FALSE), manifest(x, lp = FALSE))
    test.factor <- length(indexNew.factor)>0
  
    if(test.factor){
        if(any(endogenous(x) %in% indexOld.factor == TRUE)){
            stop("prepareData.lvm: endogenous variables must not be categorical \n",
                 "incorrect variables: ",paste(endogenous(x)[endogenous(x) %in% indexOld.factor == TRUE], collapse = " "),"\n")
        }
        ls.factor <- lapply(indexOld.factor, function(var){unique(data[[var]])})
        names(ls.factor) <- indexOld.factor
        conversion.factor <- sapply(indexNew.factor, renameFactor, ls.level = ls.factor)
    }else{
        conversion.factor <- NULL
    }
     x <- resC2D$x
     data <- resC2D$data  
  
   #### rescale data
    if(class(data)[1] != "data.frame"){data <- as.data.frame(data)}
    if(length(index.numeric)>0){
        value.center <- sapply(index.numeric, function(x){do.call(method.center,args = list(na.omit(data[[x]])))})
        value.scale <- sapply(index.numeric, function(x){do.call(method.scale,args = list(na.omit(data[[x]])))})
        data[, index.numeric] <- scale(data[, index.numeric, drop = FALSE], center = value.center, scale = value.scale)
    }
  
    #### export
    return(list(data = data,
                conversion.factor = conversion.factor,
                scale = value.scale,
                center = value.center,
                lvm = x))
}

#' @rdname initializeData
initializeData.pmultigroup <- function(x, method.center = "mean", method.scale = "sd"){

    n.model <- length(x$data)
    
    if(any(manifest(x, lp = FALSE) %in% names(x$data[[1]]) == FALSE)){
        stop("prepareData.lvm: arguments \'data\' and \'x\' are incompatible \n",
             "variables: ",paste(manifest(x)[manifest(x) %in% names(data) == FALSE], collapse = " ")," not found in \'data\' \n")
    }
  
    #### rescale data   
    for(iData in n.model){
        dataTempo <- x$data[[iData]]
        if(class(dataTempo)[1] != "data.frame"){dataTempo <- as.data.frame(dataTempo)}
        if(length(x$index.numeric)>0){
            value.center <- sapply(x$index.numeric, function(x){do.call(method.center,args = list(na.omit(dataTempo[[x]])))})
            value.scale <- sapply(x$index.numeric, function(x){do.call(method.scale,args = list(na.omit(dataTempo[[x]])))})
            dataTempo[, x$index.numeric] <- scale(dataTempo[, x$index.numeric, drop = FALSE], center = value.center, scale = value.scale)
        }
        x$data[[iData]] <- dataTempo
    }
  
    #### export
    return(list(lvm = x))
}
# }}}

# {{{ initializeOperator
#' @title Initialise a proximal operator
#' 
#' @description Generate the proximal operator for the requested penalty
#' 
#' @param lambda1 lasso penalization parameter.
#' @param lambda2 ridge penalization parameter.
#' @param lambdaG group lasso penalization parameter.
#' @param lambdaN nuclear norm penalization parameter.
#' @param index.penalty1 a list containing the position of lasso penalized parameters in the vector of parameters.
#' @param index.penalty2 a list containing the position of ridge penalized parameters in the vector of parameters.
#' @param index.penaltyG a list containing the position of group lasso parameters in the vector of parameters.
#' @param index.penaltyN a list containing the position of the parameters corresponding to each image in the vector of parameters.
#' @param index.model1 the index of the model relative to each lasso penalized parameter
#' @param index.model2 the index of the model relative to each ridge penalized parameter
#' @param index.modelG the index of the model relative to each group lasso penalized parameter
#' @param nrow a vector containing the number of rows of each image.
#' @param ncol a vector containing the number of columnss of each image.
#' @param equivariance should the lambda parameter be multiplied with the first variance parameter?
#' @param constrain.variance should the variance parameters be exponential trans
#' @param index.variance the position of the variance parameters in coef
#' @return a list containing the proximal operator and the penalty funtion
#'
#' @examples 
#' lavaPenalty:::initializeOperator(lasso = TRUE, ridge = FALSE, groupLasso = FALSE, nuclearNorm = FALSE)
#' 
initializeOperator <- function(lambda1, index.penalty1, index.model1,
                               lambda2, index.penalty2, index.model2,
                               lambdaG, index.penaltyG, index.modelG,
                               lambdaN, index.penaltyN, nrow, ncol,
                               equivariance, constrain.variance, index.variance){

  test.lasso <- length(lambda1)>0
  test.ridge <- length(lambda2)>0
  test.groupLasso <- length(lambdaG)>0
  test.nuclearNorm <- length(lambdaN)>0

  proxOperator <- list()
  objective <- list()
  
  if(test.nuclearNorm && (test.lasso || test.ridge)){
    stop("nuclear penalty in conjunction with elastic net is not yet implemented \n")
  }
  
    #### No penalty
    if(!test.lasso && !test.ridge && !test.groupLasso && !test.nuclearNorm){ 
        proxOperator$Id <- function(x, ...){x}
        objective$Id <- function(x, ...){0} 
    }else{
  
        #### Lasso penalty
        if(test.lasso){
            proxOperator$L1 <- function(x, step, lambda1, index.penalty1, ...){
                x[index.penalty1] <- mapply(proxL1, x = x[index.penalty1], step = step,
                                            lambda = lambda1)
                return(x)
            }  
    
            objective$L1 <- function(x, lambda1, index.penalty1, ...){
                return(sum(lambda1 * abs(x[index.penalty1])))
            }     
        }
  
  
        #### Ridge penalization
        if(test.ridge){
            proxOperator$L2 <- function(x, step, lambda2, index.penalty2, ...){
                x[index.penalty2] <- mapply(proxL2, x = x[index.penalty2], step = step,
                                            lambda = lambda2)
                return(x)
            }
    
            objective$L2 <- function(x, lambda2, index.penalty2, ...){
                return(sum(lambda2/2 *x[index.penalty2]^2))
            }
        }

        #### Group lasso
        if(test.groupLasso){
            # normally lambdaG has the same value for all members of a group
            # so taking the first lambda of the group is correct
            proxOperator$G <- function(x, step, lambdaG, index.penaltyG, ...){
        
                for(iGroup in 1:length(index.penaltyG)){
                    x[index.penaltyG[[iGroup]]] <- proxE2(x = x[index.penaltyG[[iGroup]]],
                                                          step = step,
                                                          lambda = lambdaG[iGroup])
                }
                return(x)
          
            }
      
            objective$G <- function(x, lambdaG, index.penaltyG, ...){
                n.group <- length(lambdaG)
                res <- lapply(1:n.group, function(g){
                    lambdaG[g]*norm(x[index.penaltyG[[g]]], type = "2")
                })
                return(sum(unlist(res)))}
        }

        #### Nuclear norm penalty
        if(test.nuclearNorm){
            proxOperator$N <-  function(x, step, lambdaN, index.penaltyN, ...){
          
          for(iNuclear in 1:length(index.penaltyN)){
              x[index.penaltyN[[iNuclear]]] <- proxNuclear(x = x[index.penaltyN[[iNuclear]]],
                                                           step = step,
                                                           lambda = lambdaN[index.penalty[[iNuclear]][1]],
                                                           nrow = nrow[iNuclear],
                                                           ncol = ncol[iNuclear])
          }
          return(x)
          
            }
            objective$N <- function(x, lambdaN, index.penaltyN, ...){
                res <- lapply(1:length(index.penaltyN), function(g){
                    lambdaN[index.penalty[[g]][1]] * sum(svd(matrix(x[index.penaltyN[[g]]], nrow = nrow[g], ncol = ncol[g]))$d)              
                })
                return(sum(unlist(res)))
            }
          
        }

    }
    
    #### compose operators
    composeOperator <- function(x, step) {

        if(!is.null(lambda1) && equivariance){
            norm1 <- x[index.variance[index.model1]]
            if(constrain.variance){
                norm1 <- exp(norm1)
            }else if(any(norm1 < 0)){
                stop("proxGrad: negative variance parameter - set constrain to TRUE in control \n")
            }
        }else{
            norm1 <- 1
        }
            
        if(!is.null(lambda2) && equivariance){
            norm2 <- x[index.variance[index.model2]]
            if(constrain.variance){
                norm2 <- exp(norm2)
            }else if(any(norm2 < 0)){
                stop("proxGrad: negative variance parameter - set constrain to TRUE in control \n")
            }
        }else{
            norm2 <- 1
        }
            
        if(!is.null(lambdaG) && equivariance){
            normG <- x[index.variance[index.modelG]]
            if(constrain.variance){
                normG <- exp(normG)
            }else if(any(normG < 0)){
                stop("proxGrad: negative variance parameter - set constrain to TRUE in control \n")
            }
        }else{
            normG <- 1
        } 

        if(constrain.variance){x[index.variance] <- exp(x[index.variance])}
        for (f in 1:length(proxOperator)) {            
            x <- proxOperator[[f]](x, step,
                lambda1 = lambda1/norm1, index.penalty1 = index.penalty1,
                lambda2 = lambda2/norm2, index.penalty2 = index.penalty2,
                lambdaG = lambdaG/normG, index.penaltyG = index.penaltyG,
                lambdaN = lambdaN,
            )            
        }
        if(constrain.variance){x[index.variance] <- log(x[index.variance])}
        return(x)
    }

    composeObjective <- function(x){
        if(equivariance){
            norm <- sum(x[index.variance])
            if(norm < 0){stop("negative variance parameter - set constrain to TRUE in control \n")}
        }else{
            norm <- 1
        }
        all.penalty <- lapply(objective, function(f){
            f(x,
              lambda1 = lambda1/norm, index.penalty1 = index.penalty1,
              lambda2 = lambda2/norm, index.penalty2 = index.penalty2,
              lambdaG = lambdaG/norm, index.penaltyG = index.penaltyG,
              lambdaN = lambdaN/norm)
        })
        return(sum(unlist(all.penalty)))
    }
    
    return(list(proxOperator = composeOperator,
                objectivePenalty = composeObjective
                ))
}

# }}}

#----------------------------------------------------------------------
### Penalty_initializer.R ends here
