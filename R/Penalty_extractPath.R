
# {{{ getPath

# {{{ documentation
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
# }}}

# {{{ getPath.plvmfit
#' @rdname getPath
#' @export
`getPath.plvmfit` <- function(x, coef, ...) {

    if(is.path(x) == FALSE){
        stop("getPath.plvmfit: no penalization path in the plvmfit object \n",
             "set argument \'regularizationPath\' to TRUE when calling estimate \n")
    }
    table.penalty <- penalty(x, type = NULL)

    
    if(!missing(coef)){
        if(identical(coef, "penalized")){
            coef <- table.penalty[penalty == "lasso",link]            
        }else if(identical(coef, "Npenalized")){
            coef <- setdiff(coef(x),table.penalty[penalty == "lasso",link])
        }
    }
    return(getPath(x$regularizationPath,
                   coef = coef,
                   V = penalty(x, type = "Vlasso")$Vlasso,
                   ...))
}
# }}}

# {{{ getPath.regPath
#' @rdname getPath
#' @export
`getPath.regPath` <- function(x, coef, V, lambda = "lambda1.abs",
                              keep.indexChange = FALSE, keep.index = TRUE, keep.optimum = TRUE, path.constrain = FALSE,
                              only.breakpoints = FALSE, increasing = TRUE, row = NULL) {

    ## extract path
    regPath <- x$path
    validNames.lambda <- c("lambda1.abs", "lambda1", "lambda2.abs", "lambda2")
    validNames.optimum <- if(!is.null(x$criterion)){c(x$criterion,"cv","optimum")}else{NULL}
    validNames.coef <- setdiff(names(regPath), c(validNames.lambda,"indexChange","index",validNames.optimum))

    # {{{ test and initialize arguments
    ## lambda
    if(missing(lambda)){
        lambda <- validNames.lambda
    }else{
        if(any(lambda %in% validNames.lambda == FALSE)){
            stop("getPath.plvmfit: invalid value for \'lambda\' \n",
                 "lambda: ",lambda,"\n",
                 "valid values: ",paste(validNames.lambda, collapse = "\" \""),"\n")
        }
    }

    ## coef
    if(missing(coef)){
        coef <- validNames.coef
    }else{
        if(any(coef %in% validNames.coef == FALSE)){
            stop("getPath.plvmfit: invalid value for \'coef\' \n",
                 "coefficient: \"",paste(coef, collapse = "\" \""),"\"\n",
                 "valid values: \"",paste(validNames.coef, collapse = "\" \""),"\"\n")
        }
    }
    
    ## only.breakpoints
    # only keep points where there is a change in the non 0 coefs (+ first and last knot))
    if(only.breakpoints){
        indexChange <- regPath[["indexChange"]]
        indexChange[is.na(indexChange)] <- -1
        test.change <- which(diff(na.omit(indexChange))!=0)+1
        row <- sort(union(c(1,NROW(regPath))
                         ,intersect(which(!is.na(regPath[["indexChange"]])), test.change)))
    }
    
    ## row    
    if(!is.null(row)){
        if(any(row %in% 1:NROW(regPath) == FALSE)){
            stop("incorrect specification of argument \'row\' \n")
        }
    }else{
        row <- 1:NROW(regPath)
    }
    # }}}
    
    # {{{ order the regularization path and convert to data.table
    seqLambda <- regPath[[lambda[1]]]
    index.order <- order(seqLambda, decreasing = 1-increasing)
    regPath <- regPath[index.order]
    seqLambda <- seqLambda[index.order]

    keep.cols <- c(if(keep.index){"index"},
                   lambda,
                   if(keep.indexChange){"indexChange"},
                   coef,
                   if(keep.optimum){validNames.optimum}
                   )
    # }}}

    # {{{ add constrain
    if(path.constrain){
        
        ## initialization
        seqKnot <- 1:NROW(regPath)
        if(getIncreasing(x)){ # start non penalized
            current.constrain <- 1:NCOL(V)
        }else{ # start fully penalized            
            current.constrain <- NULL
        }

        ## loop over the path
        ls.constrain <- NULL
        for(iKnot in seqKnot){
            ccTempo <- regPath[iKnot,indexChange]
            if(!is.na(ccTempo)){
              if(ccTempo %in% current.constrain){
                current.constrain <- setdiff(current.constrain, ccTempo)
              }else{
                current.constrain <- union(current.constrain, ccTempo)
              }
            } 
            ls.constrain <- c(ls.constrain, list(current.constrain))
        }

        regPath[, n.constrain0 := unlist(lapply(ls.constrain, length))]
        regPath[, constrain0 := ls.constrain]

        ## attribute a name to each constrain
        namesTempo <- rownames(V)
        name.penalty <- apply(V, 2, function(col){
            indexN0 <- which(col!=0)
            sign <- ifelse(col[indexN0]>0,"+","")
            sign[1] <- ""
            coef <- gsub("-1*","-",gsub("1*","",paste0(col[indexN0],"*"),fixed = TRUE), fixed = TRUE)                 
            paste(paste(sign,coef,namesTempo[indexN0],sep=""),"=0", collapse = " ",sep ="")
        })
        regPath[, constrain0.name := lapply(ls.constrain, function(c){name.penalty[c]})]
        keep.cols <- c(keep.cols, "n.constrain0","constrain0","constrain0.name")
    }
    # }}}

    ## export
    return(regPath[row,.SD,.SDcols = keep.cols])
}
# }}}

# }}}


# {{{ setPath
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
# }}}

# {{{ is.path
`is.path` <- function(x, ...) UseMethod("is.path")
`is.path.plvmfit` <- function(x) {
  return( class(x$regularizationPath) == "regPath" )
}
# }}}

# {{{ getLambda
`getLambda` <- function(x, ...) UseMethod("getLambda")
`getLambda.regPath` <- function(x, type) {

    if(type %in% c("lambda1.abs","lambda1","lambda2","lambda2.abs") == FALSE){
        stop(type," is an invalid type \n",
             "must be one of \"lambda1.abs\" \"lambda1\" \"lambda2\" \"lambda2.abs\" \n")
    }
    
  return(x$path[[type]])
}
# }}}

# {{{ getIncreasing
`getIncreasing` <- function(x, ...) UseMethod("getIncreasing")
`getIncreasing.regPath` <- function(x) {
  return(x$increasing)
}
# }}}

# {{{ getPerformance
`getPerformance` <- function(x, ...) UseMethod("getPerformance")
`getPerformance.regPath` <- function(x) {
  return(x$performance)
}
# }}}

# {{{ getOptimum
`getOptimum` <- function(x, ...) UseMethod("getOptimum")
`getOptimum.regPath` <- function(x) {
  return(x$optimum)
}
# }}}

