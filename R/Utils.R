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


formula2character <- function(x){
  
  return(gsub("[[:blank:]]","",paste(deparse(x), collapse = "+")))
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

#' @title Response variable of a formula
#' @description Return the reponse variable contained in the formula
#' 
#' @param formula a formula
#' 
#' @examples
#' \dontrun{
#' select.response(Y1~X1+X2)
#' }
select.response <- function(formula){
  return(
  setdiff(all.vars(formula),
          all.vars(delete.response(terms(formula))))
  )
}

#' @title Combine formula
#' @description Combine formula by outcome
#' 
#' @param ls.formula a list of formula
#' @param as.formula should as.formula be applied to each element of the list
#' 
#' @examples
#' \dontrun{#' 
#' combine.formula(list(Y~X1,Y~X3+X5,Y1~X2))
#' lava.options(symbol = c("~",","))
#' combine.formula(list("Y~X1","Y~X3+X5","Y1~X2"))
#' lava.options(symbol = c("<-","<->"))
#' combine.formula(list("Y<-X1","Y<-X3+X5","Y1<-X2"))
#' 
#' combine.formula(list(Y~X1,Y~X3+X1,Y1~X2))
#' combine.formula(list(Y~X1,Y~X3+X1,Y1~X2), as.unique = TRUE)
#' }
combine.formula <- function(ls.formula, as.formula = TRUE, as.unique = FALSE){
  
  if(length(ls.formula)==0){return(NULL)}
  if(class(ls.formula)=="formula"){ls.formula <- list(ls.formula)}
  
  ls.Vars <- lapply(ls.formula, initVar_link)
  
  ls.endogeneous <- unlist(lapply(ls.Vars, "[[", 1))
  ls.X <- lapply(ls.Vars, "[[", 2)
  endogenous <- unique(ls.endogeneous)
  n.endogeneous <- length(endogenous)
  
  ls.formula2 <- vector(n.endogeneous, mode = "list")
  for(iterE in 1:n.endogeneous){
    X <- unlist(ls.X[which(ls.endogeneous==endogenous[iterE])])
    if(as.unique){X <- unique(X)}
    txt <- paste(endogenous[iterE],"~",paste(X, collapse = " + "))
    ls.formula2[[iterE]] <- as.formula(txt)
  }
  
  return(ls.formula2)
}
