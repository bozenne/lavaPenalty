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
#' renameFactor <- lava.penalty:::renameFactor
#' renameFactor("x2C", ls.levels = list(x1 = 1:5, x2 = c("A","B","C")))
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


