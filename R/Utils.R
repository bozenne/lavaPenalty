# {{{ var2dummy
#' @title Convert variable names to dummy variables names
#' @description When dealing with categorical variables, the \code{estimate} function convert the categorical variables into dummy variables.
#' This function convert a set of variable names to their corresponding name in the model with dummy variables
#' @name link2dummy
#' 
#' @param x a latent variable model.
#' @param var the variable to be transformed.
#' @param data dataset according to which the model should be updated.
#' @param rm.first.factor should the first level of each categorical variable be ignored?
#' 
#' @examples
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u ~ X1+X2
#' lavaSearch2:::var2dummy(m, var = c("X1","X2"))
#' categorical(m,labels=c("M","F","MF")) <- ~X1
#' lavaSearch2:::var2dummy(m, var = c("X1","X2"))
#' categorical(m,labels=c("1","2","3")) <- ~X2
#' lavaSearch2:::var2dummy(m, var = c("X1","X2"))
#' 

`var2dummy` <-
  function(x,...) UseMethod("var2dummy")

#' @rdname link2dummy
var2dummy.list <- function(x, var, rm.first.factor = TRUE, ...){

    var <- setNames(var,var)
    ## convertion to dummy variable name for categorical variables
    factor.var <- names(x$x$attributes$labels)
    if(!is.null(var) && any(var %in% factor.var)){
        subvar <- var[var %in% factor.var]
        for(iFactor in subvar){ # iFactor <- "X1"
            newvar <- paste0(iFactor,x$x$attributes$labels[[iFactor]])
            if(rm.first.factor){newvar <- newvar[-1]}
            newvar <- setNames(newvar, rep(iFactor, length(newvar)))
            var <- c(var[names(var)!=iFactor],newvar)            
        }
    }
    return(var)
}

#' @rdname link2dummy
var2dummy.lvm <- function(x, data = NULL, ...){

    if(is.null(data)){
        data <- sim(x,1)
    }
    res <- var2dummy(lava_categorical2dummy(x, data), ...)
    return(res)
}
# }}}

# {{{ renameFactor
#' @title Link levels and coefficients
#' @description Find the categorical variable matching the name of the coefficient in lava
#'
#' @param var the name of the coefficient (binary indicator)
#' @param ls.levels a named list containing the levels for each variable
#' @param data the dataset used to find the possible categorical/factor variables
#' @param sep the character used to collapse the variable name and its levels in var
#' 
#' @examples 
#' \dontrun{
#' lavaSearch2:::renameFactor("x2C", ls.levels = list(x1 = 1:5, x2 = c("A","B","C")))
#' }
renameFactor <- function(var, ls.levels, data, sep = ""){
  
  if(!missing(data)){
    test.factor <- unlist(lapply(data, function(x){is.factor(x) + is.character(x) > 0}))
    allvars <- names(test.factor)[test.factor]
    ls.levels <- lapply(allvars, function(x){levels(data[[x]])})
  }else {
    allvars <- names(ls.levels)
  }
  test <- unlist(lapply(1:length(allvars), function(x){any(paste(allvars[x], ls.levels[[x]], sep = sep) == var)}))
  
  return(allvars[unlist(test)])
}
# }}}

# {{{ LCSseq
#' @title Common substring sequence
#' @description get the common substring sequence among a vector of strings
#' 
#' @param x a list of character
#'
#' @examples 
#' \dontrun{
#' LCSseq <- lavaSearch2:::LCSseq
#' LCSseq(list("ad","a","ad"))
#' }
LCSseq <- function(x){
  affixe <- strsplit(x[[1]], split = "")[[1]]
  
  for(iterX in 2:length(x)){
    affixe <- qualV::LCS(affixe, strsplit(x[[iterX]], split = "")[[1]])$LCS
  }
  
  return(affixe)
}
# }}}
# {{{ p2j.dgCMatrix
#' @title Extract the column index from a sparse matrix
#' @description Extract the column index from a sparse matrix
#' @name p2j
#' 
#' @param x sparse matrix
#'
#' @examples 
#' p2j <- lavaSearch2:::p2j
#' initVcoef.lvm <- lavaSearch2:::initVcoef.lvm
#' 
#' m <- lvm(Y~X1+X2+X3)
#' V <- initVcoef.lvm(m, link = coefReg(m), group = 1:length(coefReg(m)))
#' p2j(V)
#' p2j(cbind(V,0,V))
#' p2j(cbind(0,V,0,V))
#'
#' V2 <- rbind(1,V)
#' p2j(V2)
#' p2j(cbind(0,V2,0,V2))
#' 
`p2j` <-
  function(x,...) UseMethod("p2j")

#' @rdname p2j
p2j.dgCMatrix <- function(x){

    n.col <- length(x@p)-1
    nObsPerCol <-  diff(x@p)
    j <- lapply(1:n.col, function(col){rep(col, nObsPerCol[col])})
    return(unlist(j))
    
}

# }}}

