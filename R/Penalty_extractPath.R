#' @title Get the regularization path
#' @name getPath
#' @aliases getPath getPath.plvmfit
#' 
#' @description Extract the regularization path from a lvm object
#' 
#' @param x a plvmfit object or a regPath object
#' @param names names of the columns to extract
#' @param coef names of the coefficients to extract. Can also be \code{"coef0"} or \code{"coefn0"} to extract the number of 0 or non-0 coefficients
#' @param lambda names of the penalization parameter to extract.
#' @param only.breakpoints should the path be extracted only at the breakpoints
#' @param increasing should the path be extracted by increasing value of regularization parameter
#' @param order the regularization parameter to consider when ordering the points along the path
#' @param row the index at which the path should be extracted

#' @rdname getPath
#' @export
`getPath` <- function(x, ...) UseMethod("getPath")

#' @rdname getPath
#' @export
`getPath.plvmfit` <- function(x, ...) {
  
  if(isPath(x) == FALSE){
    stop("getPath.plvmfit: no penalization path in the plvmfit object \n",
         "set argument \'regularizationPath\' to TRUE when calling estimate \n")
  }
  
  getPath(x$regularizationPath, ...)
}

#' @rdname getPath
#' @export
`getPath.regPath` <- function(x, names = NULL, coefficient, lambda, only.breakpoints = FALSE, increasing = TRUE, row = NULL) {
  
  ## preparation of the request
  regPath <- x$path
  validNames <- names(regPath)
  names.penalized <- intersect(validNames, coef(x))
  name.coef <- names(regPath)[-(1:5)]
  
  if(!missing(coefficient) && !is.null(coefficient)){
    if(length(coefficient)>1){stop("getPath.plvmfit: coefficient must have length 1 \n")}
    if(coefficient %in% c("coef0", "coefn0")){only.breakpoints <- TRUE}
  }
  
  ## sort the regularization path
  regPath <- regPath[order(regPath[,getLambda(x)], decreasing = 1-increasing),,drop = FALSE]
 
  if(!is.null(names)){
    
    if(all(names %in% validNames)){
      res <- regPath[,names,drop = FALSE]
    }else{
      stop("getPath.plvmfit: ",names," contains invalid names \n",
           "invalid names: ",paste(names[names %in% validNames == FALSE], collapse = "\" \""),"\n",
           "valid names: ",paste(validNames, collapse = "\" \""),"\n")
    }
    
  } else {
    
    if(missing(lambda)){
      names.lambda <- c("lambda1.abs", "lambda1", "lambda2.abs", "lambda2")
    }else if(is.null(lambda)){
      names.lambda <- NULL
    }else{
      validValues <- c("abs", "nabs", "lambda1", "lambda2", "lambda1.abs", "lambda2.abs")
      if(any(lambda %in% validValues == FALSE)){
        stop("getPath.plvmfit: invalid value for \'lambda\' \n",
             "lambda: ",lambda,"\n",
             "valid values: ",paste(validValues, collapse = "\" \""),"\n")
      }
      
      names.lambda <- sapply(lambda, function(l){
        switch(l,
               "abs" = c("lambda1.abs", "lambda2.abs"),
               "nabs" = c("lambda1", "lambda2"),
               "lambda1" = "lambda1",
               "lambda2" = "lambda2",
               "lambda1.abs" = "lambda1.abs",
               "lambda2.abs" = "lambda2.abs"
        )
      })
    }
  
  ## extract the path
  if(only.breakpoints == 1){ # only keep points where there is a change in the non 0 coefficients (+ first and last knot)
    indexChange <- regPath$indexChange
    indexChange[is.na(indexChange)] <- -1
    test.change <- which(diff(na.omit(indexChange))!=0)+1
    index <- sort(union(c(1,NROW(regPath)),intersect(which(!is.na(regPath$indexChange)), test.change)))
    regPath <- regPath[index, ,drop = FALSE]
  }
  
  ##
  if(missing(coefficient)){
      names.coef <- setdiff(validNames, c("lambda1.abs", "lambda1", "lambda2.abs", "lambda2", "indexChange"))
    }else if(is.null(coefficient)){
      names.coef <- NULL
    }else if(coefficient %in% c("coef0","coefn0")){
     
      coefChange <- names(regPath)[regPath$indexChange+5] # 5 corresponds to the five colums "lambda1.abs", "lambda1", "lambda2.abs", "lambda2", "indexChange"
      if(getIncreasing(x)){
        current.coef <- intersect(validNames, coef(x))
        seqIterator <- 1:NROW(regPath)
      }else{
        current.coef <- NULL
        seqIterator <- NROW(regPath):1
      }
     
      ls.coef <- NULL
      for(iterPath in seqIterator){
        if(!is.na(coefChange[iterPath])){
          if(coefChange[iterPath] %in% current.coef){
            current.coef <- setdiff(current.coef, coefChange[iterPath])
          }else{
            current.coef <- union(current.coef, coefChange[iterPath])
          }
        }
        
        if(coefficient == "coefn0"){
          update <- current.coef
        }else{
          update <- setdiff(coef(x),current.coef)
        }
        
        ## add names
        if(!is.null(update)){
          attr(update, "row") <- rownames(regPath)[iterPath]
          if(length(names.lambda)>0){
            for(iterN in names.lambda){
              attr(update, iterN) <- regPath[iterPath,iterN]
            } 
          }
        }
        ls.coef <- c(ls.coef, list(update))
      }
      
      if(0 %in% regPath$lambda1.abs == FALSE){ls.coef <- rev(ls.coef)}
      return(ls.coef)
            
    } else {
      validValues <- c("all","penalized", "npenalized", "n.coef0", "n.coefn0", "coef0", "coefn0")
      if(coefficient %in% validValues == FALSE){
        stop("getPath.plvmfit: invalid value for \'coefficient\' \n",
             "coefficient: ",coefficient,"\n",
             "valid values: ",paste(validValues, collapse = "\" \""),"\n")
      }
      
      names.coef <- switch(coefficient,
                           "all" = name.coef,
                           "penalized" = names.penalized,
                           "npenalized" = setdiff(name.coef,names.penalized),
                           "n.coef0" = {regPath$coef0 <- rowSums(abs(regPath[,names.penalized, drop = FALSE] == 0)) ; "n.coef0"},
                           "n.coefn0" = {regPath$coefn0 <- rowSums(abs(regPath[,names.penalized, drop = FALSE] != 0)) ; "n.coefn0"}
      )
      
    }
    res <- regPath[,c(names.lambda, names.coef),drop = FALSE]
  }
    
  if(is.null(row)){
    return(res)
  }else if(is.numeric(row)){
    return(res[row,,drop = FALSE])
  }else if(is.character(row)){
    return(res[rownames(res) %in% row,,drop = FALSE])
  }
}

`setPath<-` <- function (x, ..., value) {
  UseMethod("setPath<-", x)
}

`setPath<-.plvmfit` <- function(x, row = NULL, names = NULL, value) {
  
  if(is.null(names)){
    names.value <- names(value)
    
    if(length(names.value)>0){
      names <- names(value)  
    }else{
      stop("setPath<-: argument names must be specified \n")
    }
  }
  
  if(is.null(row)){
    x$regularizationPath[, names] <- value
  }else{
    x$regularizationPath[row, names] <- value
  }
  
  return(x)
}

`isPath` <- function(x, ...) UseMethod("isPath")
`isPath.plvmfit` <- function(x) {
  return( class(x$regularizationPath) == "regPath" )
}

`getLambda` <- function(x, ...) UseMethod("getLambda")
`getLambda.regPath` <- function(x) {
  return(x$lambda)
}

`getIncreasing` <- function(x, ...) UseMethod("getIncreasing")
`getIncreasing.regPath` <- function(x) {
  return(x$increasing)
}

`getPerformance` <- function(x, ...) UseMethod("getPerformance")
`getPerformance.regPath` <- function(x) {
  return(x$performance)
}

`getOptimum` <- function(x, ...) UseMethod("getOptimum")
`getOptimum.regPath` <- function(x) {
  return(x$optimum)
}

`coef.regPath` <- function(x) {
  return(x$penCoef)
}
