#' @title Link levels and coefficients
#' @description Find the categorical variable matching the name of the coefficient in lava
#'
#' @param var the name of the coefficient (binary indicator)
#' @param ls.levels a named list containing the levels for each variable
#' @param data the dataset used to find the possible categorical/factor variables
#' @param sep the character used to collapse the variable name and its levels in var
#' 
#' @examples 
#' 
#' lava.penalty:::renameFactor("x2C", ls.levels = list(x1 = 1:5, x2 = c("A","B","C")))
#' 
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
#' lava.penalty:::LCSseq(list("ad","a","ad"))
#'
LCSseq <- function(x){
  affixe <- strsplit(x[[1]], split = "")[[1]]
  
  for(iterX in 2:length(x)){
    affixe <- qualV::LCS(affixe, strsplit(x[[iterX]], split = "")[[1]])$LCS
  }
  
  return(affixe)
}


#' @title Normalize var1 and var2
#' @description Convert var1 and var2 from formula or covariance to character
#' 
#' @param var1 a character indicating the endogeneous variable or a formula
#' @param var2 an optional character indicating the exogeneous variable
#' @param repVar1 should var1 be duplicated to match var2 length. Only active if format = "list".
#' @param format should the name of the variable be return (format = "list"), a vector of character formula ("txt.formula") or a list of formula ("formula")
#' 
#' @examples
#' lava.penalty:::initVar_link(var1 = a~b)
#' lava.penalty:::initVar_link(var1 = a ~ b)
#' lava.penalty:::initVar_link(var1 = a ~ b+c+d*e, format = "list")
#' lava.penalty:::initVar_link(var1 = a ~ b+c+d*e, format = "txt.formula")
#' lava.penalty:::initVar_link(var1 = a ~ b+c+d*e, format = "formula")
#' 
#' lava.penalty:::initVar_link(var1 = "a,b")
#' lava.penalty:::initVar_link(var1 = "a", var2 = "b")
#' 
#' lava.penalty:::initVar_link(var1 = Y~X1+X2)
#' lava.penalty:::initVar_link(var1 = Y~X1+X2, repVar1 = TRUE)
#' lava.penalty:::initVar_link(var1 = Y~X1+X2, format = "formula")
#' lava.penalty:::initVar_link(var1 = Y~X1+X2, format = "txt.formula")
initVar_link <- function(var1, var2, repVar1 = FALSE, format = "list"){
  
  if(missing(var2) && is.character(var1)){
    if(grepl(",",var1)==TRUE){
      var1 <- gsub(",","~", x = var1)
      sep <- ","
    }
    if(grepl("~",var1)==TRUE){
      var1 <- as.formula(var1)
      sep <- "~"
    }
  }else{
    sep <- "~"
  }
  
  if(class(var1) == "formula"){
    var2 <- all.vars(delete.response(terms(var1)))
    n.var2 <- length(var2)
    var1 <- setdiff(all.vars(var1),var2)
  }
  
  #### convert to format
  if(format == "formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){as.formula(paste(var1[i], var2[i], sep = sep))})
    
  }else if(format == "txt.formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){paste(var1[i], var2[i], sep = sep)})
    
  }else if(format == "list"){
    if(repVar1){var1 <- rep(var1, length(var2))}
    
    res <- list(var1 = var1,
                var2 = var2)
  }
 
  #### export ####
  return(res)
}

#' @title Response variable of a formula
#' @description Return the reponse variable contained in the formula
#' 
#' @param formula a formula
#' 
#' @examples
#' select.response(Y1~X1+X2)
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
#' combine.formula(list(Y~X1,Y~X3+X5,Y1~X2))
#' combine.formula(list("Y~X1","Y~X3+X5","Y1~X2"))
#' 
#' combine.formula(list(Y~X1,Y~X3+X1,Y1~X2))
#' combine.formula(list(Y~X1,Y~X3+X1,Y1~X2), as.unique = TRUE)
#' 
combine.formula <- function(ls.formula, as.formula = TRUE, as.unique = FALSE){
  
  if(length(ls.formula)==0){return(NULL)}
  
  ls.Vars <- lapply(ls.formula, function(x){
    if(as.formula){x <- as.formula(x)}
    res <- initVar_link(x)
  })
  
  
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
