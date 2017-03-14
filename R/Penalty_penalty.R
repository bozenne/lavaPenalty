
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
#' @export
penalty.lvm <- function(x,  type = NULL, nuclear = FALSE, lambdaPerCoef = FALSE, add.names = TRUE, ...){

    if(nuclear){
    if(is.null(x$penaltyNuclear$link)){return(NULL)}
    pen <- penalty(x$penaltyNuclear, type = type, ...)
    
    if(lambdaPerCoef && "lambdaN" %in% type){
      name.coef <- coef(x)
      n.coef <- length(name.coef)
      penaltyLink <- penalty(x$penalty)
      
      lambdaN <- rep(0, n.coef)
      lambdaN[name.coef %in% penaltyLink] <- penalty(x$penalty, type = "lambdaN")
      if(add.names){names(lambdaN) <- name.coef}
      if(is.list(pen)){pen$lambdaN <- lambdaN}else{pen <- lambdaN}
    }
    
  }else{
      if(is.null(x$penalty)){return(NULL)}
      pen <- penalty(x$penalty, type = type, ...)    
  }
  
  return(pen)
}
# }}}

# {{{ penalty.lvmfit
#' @rdname penaltyExtract
#' @export
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
#' @export
penalty.penaltyL12 <- function(x,
                               type = NULL,
                               no.lasso = FALSE,
                               no.ridge = FALSE,
                               no.group = FALSE,
                               index.group = NULL){ # 

    if(identical(type,"object")){

        return(x)        

    }else if(!is.null(type)){ ## extract matrix
        
        if(any(type %in% names(x) == FALSE)){
            stop("type ",paste(type[type %in% names(x) == FALSE], collapse = " ")," is not valid \n",
                 "valid types: \"",paste(names(x), collapse = "\" \""),"\" \"object\" \n")
        }

        res <- x[type]        
        return(res)
        
    }else{ ## extract table with the links
        dt.res <- NULL
        n.groups <- if(!is.null(x$Vgroup)){ NCOL(x$Vgroup) } else { 0 }

        ## lasso penalty
        if(!is.null(x$Vlasso) && no.lasso == FALSE){
            dt.lasso <- data.table::data.table(link = x$Vlasso@Dimnames[[1]][x$Vlasso@i+1],
                                               group = p2j(x$Vlasso),
                                               coef = x$Vlasso@x,
                                               penalty = "lasso")
            dt.res <- rbind(dt.res, dt.lasso)
        }
        ## ridge penalty
        if(!is.null(x$Vridge) && no.ridge == FALSE){
            dt.ridge <- data.table::data.table(link = x$Vridge@Dimnames[[1]][x$Vridge@i+1],
                                               group = p2j(x$Vridge),
                                               coef = x$Vridge@x,
                                               penalty = "ridge")
            dt.res <- rbind(dt.res, dt.ridge)
        }
        ## group penalty
        if(n.groups > 0 && no.group == FALSE){

            dt.group <- data.table::data.table(link = x$Vgroup@Dimnames[[1]][x$Vgroup@i+1],
                                               group = p2j(x$Vgroup),
                                               coef = x$Vgroup@x,
                                               penalty = "group lasso")
            if(!is.null(index.group)){
                if(any(group %in% 1:n.groups == FALSE) && any(group %in% dt.group[["group"]] == FALSE)){
                    stop("unknown groups ",paste(group[group %in% 1:n.groups == FALSE], collapse = " ")," \n")
                }
                dt.res <- rbind(dt.res, dt.group[group %in% index.group])
            }else{
                dt.res <- rbind(dt.res, dt.group)
            }
            
        }

        if(!is.null(dt.res)){
            resInit <- lava.reduce::initVar_links(dt.res$link, format = "list")
            dt.res[, endogenous := resInit$var1]
            dt.res[, exogenous := resInit$var2]
        }

        ## export
        return(dt.res)
    }
    

}
# }}}

# {{{ penalty.penaltyNuclear
#' @rdname penaltyExtract
#' @export
penalty.penaltyNuclear <- function(x, type, group = NULL, keep.list = FALSE){
  
     if(identical(type,"object")){
        return(x)
    }else if(is.null(type)){
        res <- x
        type <- names(x)
    }else if(any(type %in% names(x) == FALSE)){
        stop("type ",paste(type[type %in% names(x) == FALSE], collapse = " ")," is not valid \n",
             "valid types: \"",paste(names(x), collapse = "\" \""),"\" \"object\" \n")
    }else{
        res <- x[type]
    }

    return(res)
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
#' pm
#' 
#' @export
`penalty<-` <-
  function(x,...) UseMethod("penalty<-", x)

# }}}

# {{{ penalty<-.plvm
#' @rdname penaltyUpdate
#' @export
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
#' @export
`penalty<-.penaltyL12` <- function(x, type, lambda, add, value){

    validTypes <- names(x)
    
    if(is.null(type)){
    
    x <- value
    
    }else if(length(type)==1){

        ## check type
        if(type %in% names(x) == FALSE){
            stop("type ",type," is not valid \n",
                 "valid types: \"",paste(names(x), collapse = "\" \""),"\" \n")
        }

        ## combine with existing values        
        if(add==TRUE && length(x[[type]])>0){
            if(is.matrix(value)){
                value <- cbind(x[[type]],value)
            }else{
                value <- c(x[[type]],value)
            }
        }

        ## check matching V-lambda
        type2 <- switch(type,
                        "lambda1" = "Vlasso",
                        "lambda2" = "Vridge",
                        "lambdaG" = "Vgroup",
                        "Vlasso" = "lambda1",
                        "Vridge" = "lambda2",
                        "Vgroup" = "lambdaG"
                        )

        if(type %in% paste0("lambda",c("1","2","G"))){
            n.link <- length(x[[type2]]@x)
            if(length(value)==1){                
                value <- rep(value, n.link)
            }else if(length(value) != n.link){
                stop("\'",type,"\' must have length ",n.link,"\n",
                     "current length: ",length(value),"\n")
            }
        }else{
            if(unique(length(x[[type2]])) == 1){
                n.link <- length(x[[type]]@x)
                x[[type2]] <- rep(x[[type2]][1], n.link)
            }else{
                x[[type2]] <- numeric(0)
            }
        }

        ## affect value
        x[[type]] <- value
        
    }else{
    
        stop("argument \'type\' must have length 0 or 1 \n")
    
    }
  
  
  return(x)
}
# }}}

# {{{ penalty<-.penaltyNuclear
#' @rdname penaltyUpdate
#' @export
`penalty<-.penaltyNuclear` <- function(x, type = "link", value){
  browser()
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
