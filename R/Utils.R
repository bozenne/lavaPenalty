#' @title Find the categorical variable matching the name of the coefficient in lava
#' @examples 
#' renameFactor("x2C", ls.levels = list(x1 = 1:5, x2 = c("A","B","C")))
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
  return(deparse(x))
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