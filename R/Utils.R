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
