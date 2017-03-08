
# {{{ penalty

# {{{ doc
#' @title Extract the penalty from a lvm object
#' @name penaltyExtract
#' @description Extract the penalty from a lvm object
#'
#' @param x a lvm model
#' @param type the information about the penalty to be returned.
#' @param group the group of penalty parameter to be returned
#' @param nuclear should informations be extracted for the nuclear penalty
#' @param keep.list should the format of the results be always a list
#' @param lambdaPerCoef when a regularization parameter is extracted should it expanded to match the parameters of the LVM
#' @param add.names when a regularization parameter is extracted should the name of the parameters be added in the output. Only active if lambdaPerCoef is TRUE.
#' 
#' @examples 
#' set.seed(10)
#' n <- 300
#' formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:12), collapse = "+")))
#' mSim <- lvm(formula.lvm)
#' df.data <- sim(mSim,n)
#' 
#' pm <- penalize(mSim, c("Y~X1","Y~X4", "Y~X10"))
#' penalty(pm, type = "link")
#' penalty(pm, type = c("link","group"))
#' 
#' pm <- penalize(pm, value = paste0("Y~X",5:6), group = 1, add = TRUE)
#' penalty(pm, type = "link") # all parameters
#' penalty(pm, type = "link", group = 0) # individually penalized parameters
#' penalty(pm, type = "link", group = 1) # group 1 penalization
#' penalty(pm, type = "V")
#' penalty(pm, type = "V", group = 0)
#' 
#' penalty(pm, type = "lambda1")
#' penalty(pm, type = "lambda1", lambdaPerCoef = TRUE)
#' penalty(pm, type = "lambda1", lambdaPerCoef = TRUE, add.names = FALSE)
#' penalty(pm, type = c("lambda1","lambda2"), lambdaPerCoef = TRUE)
#' 
#' penalty(pm, type = "link") <- c("Y~X1" , "Y~X4")
#' 
#' @export
`penalty` <-
  function(x,...) UseMethod("penalty")
# }}}

# {{{ penalty.lvm
#' @rdname penaltyExtract
penalty.lvm <- function(x,  type = "link", nuclear = FALSE, lambdaPerCoef = FALSE, add.names = TRUE, ...){

  if(nuclear){
    if(is.null(x$penaltyNuclear)){return(NULL)}
    
    pen <- penalty(x$penaltyNuclear, type = type, ...)
    
    if(lambdaPerCoef && "lambdaN" %in% type){
      name.coef <- coef(x)
      n.coef <- length(name.coef)
      penaltyLink <- penalty(x$penalty, type = "link")
      
      lambdaN <- rep(0, n.coef)
      lambdaN[name.coef %in% penaltyLink] <- penalty(x$penalty, type = "lambdaN")
      if(add.names){names(lambdaN) <- name.coef}
      if(is.list(pen)){pen$lambdaN <- lambdaN}else{pen <- lambdaN}
    }
    
  }else{
    if(is.null(x$penalty)){return(NULL)}

    pen <- penalty(x$penalty, type = type, ...)

    if(lambdaPerCoef && any(c("lambda1","lambda2") %in% type) ){
      name.coef <- coef(x)
      n.coef <- length(name.coef)
      penaltyLink <- penalty(x$penalty, type = "link")
      
      if("lambda1" %in% type){
        lambda1 <- rep(0, n.coef)
        lambda1[name.coef %in% penaltyLink] <- penalty(x$penalty, type = "lambda1")
        if(add.names){names(lambda1) <- name.coef}
        if(is.list(pen)){pen$lambda1 <- lambda1}else{pen <- lambda1}
      }
      if("lambda2" %in% type){
        lambda2 <- rep(0, n.coef)
        lambda2[name.coef %in% penaltyLink] <- penalty(x$penalty, type = "lambda2")
        if(add.names){names(lambda2) <- name.coef}
        if(is.list(pen)){pen$lambda2 <- lambda2}else{pen <- lambda2}
      }
    }
    
  }
  
  return(pen)
}
# }}}

# {{{ penalty.lvmfit
#' @rdname penaltyExtract
penalty.lvmfit <- function(x, type, ...){
  
  type.origin <- type
  if("lambda1.abs" %in% type){
  type[type == "lambda1.abs"] <- "lambda1"
  }
  if("lambda2.abs" %in% type){
    type[type == "lambda2.abs"] <- "lambda2"
  }
  res <- penalty.lvm(x, ...)
  
  if("lambda1.abs" %in% type){
    browser()
  }
  if("lambda2.abs" %in% type){
    browser()
  }
  return(res)
}
# }}}
# {{{ penalty.penaltyL12
#' @rdname penaltyExtract
penalty.penaltyL12 <- function(x,
                               type,                               
                               group = NULL,
                               no.elasticNet = FALSE,
                               no.group = FALSE,
                               keep.list = FALSE){ # 

    ## check and initialize arguments    
    penaltyType <- names(x)#c("lambda1","lambda2","adaptive","objectivePenalty","VelasticNet","Vgroup")  
    valideType <- c("link", "group", penaltyType)
    x$link <- NULL ; x$group <- NULL
    
    if(is.null(type)){
        type <- names(x)
    }else{
        if(any(type %in% valideType == FALSE)){
            stop("type ",paste(type[type %in% valideType == FALSE], collapse = " ")," is not valid \n",
                 "valid types: \"",paste(valideType, collapse = "\" \""),"\" \n")
        }
    }
    n.groups <- if(!is.null(x$Vgroup)){ NCOL(x$Vgroup) } else { 0 }
    if(is.numeric(group) && any(group %in% 1:n.groups == FALSE)){
        stop("unknown groups ",paste(group[group %in% 1:n.groups == FALSE], collapse = " ")," \n",
             "valid groups: ",paste(1:n.groups,collaspe = " "),"\n")
    }
    
    ## extract penalty
    if(no.elasticNet == FALSE){
        if(is.null(x$VelasticNet)){return(NULL)}
        if(any(duplicated(x$VelasticNet@i))){
            warning("multiple coefficients penalized in a single penalty \n",
                    "vector of penalised coefficients will be extracted regardless their potential interaction with other coefficients \n")
        }
        if("link" %in% type){
            x$link <- x$VelasticNet@Dimnames[[1]][unique(x$VelasticNet@i)+1]
        }
        if("group" %in% type){
            x$group <- rep(NA,length(x$link))
        }
    }else{
        if("link" %in% type){ x$link <- NULL }
        if("group" %in% type){ x$group <- NULL }
    }

    ## restrict to group of penalty    
    if(is.null(group)){
        if(n.groups == 0){ group <- FALSE } else { group <- 1:n.groups }
    }

    if(no.group == FALSE && n.groups > 0){

        if(any(group %in% 1:NCOL(x$Vgroup))){
            if("link" %in% type){
                x$link <- c(x$link,
                            x$Vgroup@Dimnames[[1]][x$Vgroup[,group,drop=FALSE]@i+1]
                            )
            }
            if("group" %in% type){
                x$group <- c(x$group,
                             unlist(sapply(group, function(g){
                                 rep(g, times = diff(x$Vgroup@p)[g])
                             }))
                             )
            }
            if("Vgroup" %in% type){
                x$Vgroup <- x$Vgroup[,group,drop=FALSE]
            }
        }        
    }

    ## export
    res <- x[type]
    if(identical(type,penaltyType)){ # keep class
        class(res) <- class(x)
    }else if((keep.list == FALSE) && length(type)==1){ # convert to vector
        res <- res[[1]]
    }

    return(res)  
}
# }}}

# {{{ penalty.penaltyNuclear
#' @rdname penaltyExtract
penalty.penaltyNuclear <- function(x, type, group = NULL, keep.list = FALSE){
  
  valideType <-  names(x)
  
  ## extract penalty
  if(is.null(type)){
    res <- x
    type <- names(x)
  }else if(any(type %in% valideType == FALSE)){
    stop("type ",paste(type[type %in% valideType == FALSE], collapse = " ")," is not valid \n",
         "valid types: \"",paste(valideType, collapse = "\" \""),"\" \n")
  }else{
    res <- x[type]
  }
  
  ## export
  if((keep.list == FALSE) && length(type)==1){
    return(res[[1]])
  }else{
    return(res)
  }
}
# }}}

# }}}

# {{{ penatly<-
# {{{ doc
#' @title Update the penalty term
#' @name penaltyUpdate
#' @description Update the penalty term of a lvm.penalty object
#'
#' @param x a lvm.penalty object
#' @param type the type of penalty information to be updated
#' @param nuclear should informations be extracted for the nuclear penalty
#' @param value the value to be attributed
#' 
#' @examples 
#' set.seed(10)
#' n <- 300
#' formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:12), collapse = "+")))
#' mSim <- lvm(formula.lvm)
#' df.data <- sim(mSim,n)
#' 
#' pm <- penalize(mSim, c("Y~X1","Y~X4", "Y~X10"))
#' 
#' penalty(pm, type = "link") <- c("Y~X5")
#' 
#' @export
`penalty<-` <-
  function(x,...) UseMethod("penalty<-", x)
# }}}

# {{{ penalty<-.plvm
#' @rdname penaltyUpdate
`penalty<-.plvm` <- function(x, nuclear = FALSE, value, ...){
  
  if(nuclear){
    pen <- x$penaltyNuclear
  }else{
    pen <- x$penalty
  }
  
  penalty(pen, ...) <- value
  
  if(nuclear){
    x$penaltyNuclear <- pen
  }else{
    x$penalty <- pen
  }
  
  return(x)
}
# }}}
# {{{ penalty<-.penaltyL12
#' @rdname penaltyUpdate
`penalty<-.penaltyL12` <- function(x, type = "link", value){
  
  validTypes <- names(x)
  
  if(is.null(type)){
    
    x <- value
    
  }else if(length(type)==1){
    
    if(type %in% validTypes == FALSE){
      stop("type ",type," is not valid \n",
           "valid types: \"",paste(validTypes, collapse = "\" \""),"\" \n")
    }
    x[[type]] <- value
    
  }else{
    
    if(any(type %in% validTypes == FALSE)){
      stop("type \"",paste(type[type %in% validTypes == FALSE], collapse =  "\" \""),"\" is not valid \n",
           "valid types: \"",paste(validTypes, collapse = "\" \""),"\" \n")
    }
    x[type] <- value
    
  }
  
  
  return(x)
}
# }}}

# {{{ penalty<-.penaltyNuclear
#' @rdname penaltyUpdate
`penalty<-.penaltyNuclear` <- function(x, type = "link", value){
  
  validTypes <- names(x)
  
  if(is.null(type)){
    
    x <- value
    
  }else if(length(type)==1){
    
    if(type %in% validTypes == FALSE){
      stop("type ",type," is not valid \n",
           "valid types: \"",paste(validTypes, collapse = "\" \""),"\" \n")
    }
    x[[type]] <- value
    
  }else{
    
    if(any(type %in% validTypes == FALSE)){
      stop("type \"",paste(type[type %in% validTypes == FALSE], collapse =  "\" \""),"\" is not valid \n",
           "valid types: \"",paste(validTypes, collapse = "\" \""),"\" \n")
    }
    x[type] <- value
    
  }
  
  return(x)
}
# }}}

# }}}
