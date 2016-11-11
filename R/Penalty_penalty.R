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

#' @rdname penaltyExtract
penalty.lvm <- function(x,  type = "link", nuclear = FALSE, lambdaPerCoef = FALSE, add.names = TRUE, ...){
  
  if(is.null(x$penalty) && is.null(x$penaltyNuclear)){ ## no penalisation
    pen <- NULL
  }else if(nuclear){
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

#' @rdname penaltyExtract
penalty.penaltyL12 <- function(x, type, group = NULL, keep.list = FALSE){

  valideType <- c("link","lambda1","lambda2",
                  "adaptive","objectivePenalty","proxOperator","Vlasso","Vridge","Vgroup")
  
  ## extract penalty
  if(is.null(type)){
    res <- x
    type <- names(x)
  }else if(any(type %in% valideType == FALSE)){
    stop("type ",paste(type[type %in% valideType == FALSE], collapse = " ")," is not valid \n",
         "valid types: \"",paste(valideType, collapse = "\" \""),"\" \n")
  }else{
    if("link" %in% type){
      if(is.null(x$V)){return(NULL)}
      yNames <- x$V@Dimnames[[1]][unique(x$V@i)+1]
      xNames <- unlist(mapply(rep, x  = x$V@Dimnames[[2]], times = diff(x$V@p)))
      x$link <- paste(yNames, tapply(xNames, x$V@i, paste, collapse = "+"), sep = lava.options()$symbols[1]) 
    }
    res <- x[type]
  }
  
  ## restrict to group of penalty
  if(!is.null(group)){
    
    indexGroup <- which(floor(x[["group"]]) == group)
    
    if("link" %in% type){
      res[["link"]] <- res[["link"]][indexGroup, drop = FALSE]
    }
    if("group" %in% type){
      res[["group"]] <- res[["group"]][indexGroup, drop = FALSE]
    }
    if("V" %in% type){
      nameGroup <- x[["link"]][indexGroup]
      res[["V"]] <- res[["V"]][nameGroup, nameGroup, drop = FALSE]
    }
    
  }
  
  ## export
  if((keep.list == FALSE) && length(type)==1){
    return(res[[1]])
  }else{
    return(res)
  }
  
}

#' @rdname penaltyExtract
penalty.penaltyNuclear <- function(x, type, group = NULL, keep.list = FALSE){
  
  valideType <- c("link","lambdaN","name.Y","name.X","nrow","ncol")
  
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

#' @rdname penaltyUpdate
`penalty<-.penaltyL12` <- function(x, type = "link", value){
  
  validTypes <- c("link","lambda1","lambda2","adaptive","objectivePenalty","proxOperator","Vlasso","Vridge","Vgroup")
  
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

#' @rdname penaltyUpdate
`penalty<-.penaltyNuclear` <- function(x, type = "link", value){
  
  validTypes <- c("link","lambdaN","name.Y","name.X","nrow","ncol")
  
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

