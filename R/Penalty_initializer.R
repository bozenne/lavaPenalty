### Penalty_initializer.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar  6 2017 (11:51) 
## Version: 
## last-updated: mar  7 2017 (15:38) 
##           By: Brice Ozenne
##     Update #: 44
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ doc

#' @title Intialization of the parameter of the plvm
#' @name initializer
#' @description Compute the coefficients of the non and the completely regularized latent variable model.
#' Also prepare the penalty.
#' 
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param name.coef the name of the parameters
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- y1~x1+x2+x3+x4
#'
#' #### elastic net
#' e <- estimate(m, sim(m, 1e2))
#' pm <- penalize(m, lambda1 = 2, lambda2 = 1.5)#' 
#' pen12 <- initializer(pm$penalty, name.coef = names(coef(e)))
#' pen12
#'
#' i <- initializer(pm, data = sim(m, 1e2),
#'                  regularizationPath = FALSE, constrain.variance = FALSE)
#' 
#' #### group lasso
#' categorical(m, labels = c("A","B","C")) <- "x1"
#' categorical(m, labels = c("A","B","C")) <- "x2"
#'
#' e <- estimate(m, sim(m, 1e2))
#' pm <- penalize(m, lambdaG = 5)
#' pen12 <- initializer(pm$penalty, name.coef = names(coef(e)))
#' pen12
#' estimate(pm, data = sim(pm, 1e2))
#'
#' i <- initializer(pm, data = sim(m, 1e2),
#'                  regularizationPath = FALSE, constrain.variance = FALSE)
#' 
#' # nuclear norm
#' m <- lvm()
#' coords <- expand.grid(x = 1:10, y = 1:10)
#' m <- regression(m, y = "y1", x = paste0("z",1:100))
#'
#' mNuclear <- lvm(y1 ~ x1 + x2)
#' penalizeNuclear(mNuclear, coords = coords, lambdaN = 10) <- coefReg(m, value = TRUE)
#' mNuclear
#'
#' penN <- initializer(mNuclear$penaltyNuclear, name.coef = coef(mNuclear))
#' penN
#' 
`initializer` <-
  function(x,...) UseMethod("initializer")

# }}}

# {{{ initializer.plvm

#' @rdname initializer
initializer.plvm <- function(x, data, regularizationPath,
                             constrain.variance, name.variance,
                             ...){

    link.penalty <- c(penalty(x, type = "link"),
                      penalty(x, type = "link", nuclear = TRUE))
    test.regularizationPath <- !is.null(regularizationPath)
    class(x) <- setdiff(class(x),"plvm")
    data <- data[,manifest(x, lp = FALSE, xlp = TRUE), drop = FALSE]
  
    names.coef <- coef(x)
    # names.coef  <- coef(lava::fixsome(x, measurement.fix = TRUE, n = NROW(data), debug = FALSE))
    n.coef <- length(names.coef)
    n.data <- NROW(data)

    # {{{ compute start - low dimensional (keep all links unpenalized)
    if(test.regularizationPath && n.data > n.coef){
        suppressWarnings(
            initLVM <- tryCatch(estimate(x = x, data = data, quick = TRUE, ...),
                                error = function(e){e}
                                )
        )
        if(all("error" %in% class(initLVM) == FALSE) ){ # check convergence
            start_lambda0 <- initLVM                 
        }else{ 
            start_lambda0 <- NULL
        }
    }else{ 
        start_lambda0 <- NULL
    }

    # }}}

    # {{{ compute start - high dimensional (remove penalized links, all lambda are infinity)

    start_lambdaMax <- setNames(rep(0,n.coef),names.coef)

    ## simplify the LVM
    x0 <- x
    
    ## remove penalized regression links
    # lasso/ridge
    reg.penalty <- intersect(coefReg(x, value = TRUE),
                             penalty(x, type = "link", no.group = TRUE))
    if(!is.null(reg.penalty)){
        ls.formula <- lava.reduce::combine.formula(reg.penalty)
        for(iter_link in ls.formula){
            x0 <- rmLink(x0, iter_link, simplify = TRUE)
        }
    }
    # group lasso
    reg.penalty <- penalty(x, type = "link", no.elasticNet = TRUE)
    Mcorrespondance <- initGroup.lvm(x, link = coefReg(x, value = TRUE), data = data)$Mcat
    
    if(!is.null(reg.penalty)){
        # in addition to the penalties remove the categorical regression link that factorize the penalties
        add.reg <- Mcorrespondance[rownames(Mcorrespondance) %in% reg.penalty,][["originalLink"]]
        ls.formula <- lava.reduce::combine.formula(c(reg.penalty, unique(add.reg)))
        for(iter_link in ls.formula){
            x0 <- rmLink(x0, iter_link, simplify = TRUE)
        }
    }
    
    # nuclearNorm
    reg.penalty <- intersect(coefReg(x, value = TRUE),
                             penalty(x, type = "link", nuclear= TRUE))
    if(!is.null(reg.penalty)){
        ls.formula <- lava.reduce::combine.formula(reg.penalty)
        for(iter_link in ls.formula){
            x0 <- rmLink(x0, iter_link, simplify = TRUE)
        }
    }
    
    # remove penalised LV [TODO]
    
    ## estimate the model
    suppressWarnings(
        newCoef <- estimate(x = x0, data = data, quick = TRUE, ...)
    )
    newCoef <- newCoef[names(newCoef) %in% names(start_lambdaMax)] # only keep relevant parameters
    start_lambdaMax[names(newCoef)] <- newCoef
  
    # }}}
 
    # export
    return(list(lambda0 = start_lambda0,
                lambdaMax = start_lambdaMax))
}

# }}}

# {{{ initializer.penaltyL12
#' @rdname initializer
initializer.penaltyL12 <- function(x, name.coef){

    n.coef <- length(name.coef)

    ## check
    name.allPenalty <- penalty(x, type = "link", no.group = FALSE, no.elasticNet = FALSE)
    if(any(name.allPenalty %in% name.coef == FALSE)){
        stop("initPenalty: some penalty will not be applied because the corresponding parameter is used as a reference \n",
             "non-applied penalty: ",paste(setdiff(name.allPenalty, name.coef), collapse = " "),"\n")
    }


    #### elastic net penalty
    lambda1 <- setNames(rep(0, n.coef), name.coef)
    lambda2 <- setNames(rep(0, n.coef), name.coef)

    name.elasticNet <- penalty(x, type = "link", no.group = TRUE, no.elasticNet = FALSE)
    if(length(name.elasticNet)>0){
            
        lambda1[name.elasticNet] <- penalty(x, type = "lambda1")        
        lambda2[name.elasticNet] <- penalty(x, type = "lambda2")

    }

    #### group lasso penalty
    lambdaG <- setNames(rep(0, n.coef), name.coef)
    index.penaltyG <- list()

    name.groupLasso <- penalty(x, type = "link", no.group = FALSE, no.elasticNet = TRUE)
    if(length(name.groupLasso)>0){
        lambdaG[name.groupLasso] <- penalty(x, type = "lambdaG")

        index.tempo <- setNames(penalty(x, type = "group"), name.groupLasso)
        index.penaltyG <- tapply(index.tempo,index.tempo, function(x){
            list(setNames(match(names(x), name.coef), names(x)))
        })

    }

    return(list(lambda1 = lambda1,
                lambda2 = lambda2,
                lambdaG = lambdaG,
                index.penaltyG = index.penaltyG))
}
# }}}

# {{{ initializer.penaltyNuclear
#' @rdname initializer
initializer.penaltyNuclear <- function(x, name.coef){
    n.coef <- length(name.coef)
    name.allPenalty <- penalty(x, type = "link")
    lambdaN <- penalty(x, type = "lambdaN")
    if(length(name.allPenalty)==0){
        return(list(lambdaN = 0, ncol = NULL, nrow = NULL, index.penaltyN = NULL))
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



# {{{ initializeOperator
#' @title Initialise a proximal operator
#' 
#' @description Generate the proximal operator for the requested penalty
#' 
#' @param lambda1 lasso penalization parameter.
#' @param lambda2 ridge penalization parameter.
#' @param lambdaG group lasso penalization parameter.
#' @param index.penaltyG a list containing the position of each group of parameters in the vector of parameters.
#' @param lambdaN nuclear norm penalization parameter.
#' @param index.penaltyN a list containing the position of the parameters corresponding to each image in the vector of parameters.
#' @param nrow a vector containing the number of rows of each image.
#' @param ncol a vector containing the number of columnss of each image.
#' @param constrain.lambda should the lambda parameter be normalized by the total variance?
#' @param constrain.variance should the variance parameters be exponential trans
#' @param index.variance
#' @return a list containing the proximal operator and the penalty funtion
#'
#' @examples 
#' lava.penalty:::initializeOperator(lasso = TRUE, ridge = FALSE, groupLasso = FALSE, nuclearNorm = FALSE)
#' 
initializeOperator <- function(lambda1,
                               lambda2,
                               lambdaG, index.penaltyG,
                               lambdaN, index.penaltyN, nrow, ncol,
                               constrain.lambda, constrain.variance, index.variance){

  test.lasso <- any(lambda1 > 0)
  test.ridge <- any(lambda2 > 0)
  test.groupLasso <- any(lambdaG > 0)
  test.nuclearNorm <- any(lambdaN > 0)

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
            proxOperator$L1 <- function(x, step, lambda1, ...){
                x[lambda1>0] <- mapply(proxL1, x = x[lambda1>0], step = step,
                                       lambda = lambda1[lambda1>0])
                return(x)
            }  
    
            objective$L1 <- function(x, lambda1, ...){
                return(sum(lambda1 * abs(x)))
            }     
        }
  
  
        #### Ridge penalization
        if(test.ridge){
            proxOperator$L2 <- function(x, step, lambda2, ...){
                x[lambda2>0] <- mapply(proxL2, x = x[lambda2>0], step = step,
                                       lambda = lambda2[lambda2>0])
                return(x)
            }
    
            objective$L2 <- function(x, lambda2, ...){
                return(sum(lambda2/2 *x^2))
            }
        }

        #### Group lasso
        if(test.groupLasso){
            # normally lambdaG has the same value for all members of a group
            # so taking the first lambda of the group is correct
            proxOperator$G <- function(x, step, lambdaG, ...){
        
                for(iGroup in 1:length(index.penaltyG)){
                    x[index.penaltyG[[iGroup]]] <- proxE2(x = x[index.penaltyG[[iGroup]]],
                                                          step = step,
                                                          lambda = lambdaG[index.penaltyG[iGroup][1]])
                }
                return(x)
          
            }
      
            objective$G <- function(x, lambdaG, ...){
                res <- lapply(index.penaltyG, function(g){lambdaG[g[1]]*norm(x[g], type = "2")})
                return(sum(unlist(res)))}
        }

        #### Nuclear norm penalty
        if(test.nuclearNorm){
            proxOperator$N <-  function(x, step, lambdaN, ...){
          
          for(iNuclear in 1:length(index.penaltyN)){
              x[index.penaltyN[[iNuclear]]] <- proxNuclear(x = x[index.penaltyN[[iNuclear]]],
                                                           step = step,
                                                           lambda = lambdaN[index.penalty[[iNuclear]][1]],
                                                           nrow = nrow[iNuclear],
                                                           ncol = ncol[iNuclear])
          }
          return(x)
          
            }
            objective$N <- function(x, lambdaN, ...){
                res <- lapply(1:length(index.penaltyN), function(g){
                    lambdaN[index.penalty[[g]][1]] * sum(svd(matrix(x[index.penaltyN[[g]]], nrow = nrow[g], ncol = ncol[g]))$d)              
                })
                return(sum(unlist(res)))
            }
          
        }

    }
    
    #### compose operators
    composeOperator <- function(x, step) {

        if(constrain.lambda){
            if(constrain.variance){
                norm <- sum(exp(x[index.variance]))
            }else{
                norm <- sum(x[index.variance])
                if(norm < 0){stop("proxGrad: negative variance parameter - set constrain to TRUE in control \n")}
            }          
        }else{
            norm <- 1
        }

        if(constrain.variance){x[index.variance] <- exp(x[index.variance])}
        for (f in proxOperator) {
             x <- f(x, step,
                    lambda1 = lambda1/norm,
                    lambda2 = lambda2/norm,
                    lambdaG = lambdaG/norm,
                    lambdaN = lambdaN/norm)
        }
        if(constrain.variance){x[index.variance] <- log(x[index.variance])}
        return(x)
    }

    composeObjective <- function(x){

        if(constrain.lambda){
            norm <- sum(x[index.variance])
            if(norm < 0){stop("negative variance parameter - set constrain to TRUE in control \n")}
        }else{
            norm <- 1
        }

        all.penalty <- lapply(objective, function(f){
            f(x,
              lambda1 = lambda1/norm,
              lambda2 = lambda2/norm,
              lambdaG = lambdaG/norm,
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
