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
#' lava.penalty:::renameFactor("x2C", ls.levels = list(x1 = 1:5, x2 = c("A","B","C")))
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
#' LCSseq <- lava.penalty:::LCSseq
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
#' p2j <- lava.penalty:::p2j
#' initVcoef.lvm <- lava.penalty:::initVcoef.lvm
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

