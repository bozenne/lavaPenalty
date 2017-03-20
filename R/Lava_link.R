# {{{ findNewLink
#' @title Find the new possible links between variables (copied from lava::modelsearch)
#' @name findNewLink
#' 
#' @param x a lvm model
#' @param rm.latent ignore links between latent variables
#' @param rm.endo ignore links between endogenous variables
#' @param output return the names of the variables to link ("names") or their position ("index")
#' 
#' @examples 
#' \dontrun{
#' findNewLink <- lava.penalty:::findNewLink
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' findNewLink(m, rm.endo = FALSE)
#' findNewLink(m, rm.endo = TRUE)
#' 
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' 
#' findNewLink(m, rm.endo = FALSE)
#' findNewLink(m, rm.endo = TRUE)
#' findNewLink(m, rm.endo = TRUE, output = "index")
#' }
`findNewLink` <-
  function(x,...) UseMethod("findNewLink")

#' @rdname findNewLink
findNewLink.lvm <- function(x, rm.latent = FALSE, rm.endo = FALSE, output = "names"){
  
  match.arg(output, choices = c("names","index"))
  
  AP <- with(lava::index(x), A + t(A) + P)
  
  restricted <- c()
  for (i in seq_len(ncol(AP) - 1)){
    for (j in seq(i + 1, nrow(AP))){
      test.latent <- (rownames(AP)[i] %in% latent(x)) + (colnames(AP)[j] %in% latent(x))
      test.endo2 <- (rownames(AP)[i] %in% endogenous(x)) + (colnames(AP)[j] %in% endogenous(x))
      
      if (AP[j, i] == 0 && (rm.latent == FALSE || test.latent!=2) && (rm.endo == FALSE || test.endo2!=2)){
        restricted <- rbind(restricted, c(i, j))
      }
    }
  }
  
  
  if (is.null(restricted)){
    return(NULL)
  }
  
  if(output == "names"){
    names.restricted <- cbind(rownames(AP)[restricted[,1]],
                              colnames(AP)[restricted[,2]])
    
    return(names.restricted)
  }else{
      return(restricted)
  }
  
}
# }}}


# {{{ addLink
#' @title Add a new link between two variables in a lvm
#' @rdname addLink
#' @description Generic interface to add links to a lvm.
#' 
#' @param x a lvm model
#' @param var1 the first variable (character) or a formula describing the link to be added to the lvm
#' @param var2 the second variable (character). Only used if var1 is a character.
#' @param covariance does the link is bidirectional. Ignored if one of the variable is exogenous.
#' @param warnings should a warning be displayed when no link is added
#' 
#' @examples 
#' \dontrun{
#' addLink <- lava.penalty:::addLink
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' m2 <- m
#' 
#' addLink(m, x1 ~ y1)
#' addLink(m, y1 ~ x1)
#' coef(addLink(m, y1 ~ y2, covariance = TRUE))
#' 
#' addLink(m2, "x1", "y1")
#' addLink(m2, "y1", "x1")
#' coef(addLink(m, "y1", "y2", covariance = TRUE))
#' }
`addLink` <-
    function(x,...) UseMethod("addLink")

#' @rdname addLink
addLink.lvm <- function(x, var1, var2, covariance, allVars = vars(x), warnings = FALSE){
  
  res <- initVar_link(var1, var2, format = "list")
  var1 <- res$var1
  var2 <- res$var2
  
  if(var1 %in% allVars == FALSE){
    if(warnings){
      warning("addLink.lvm: var1 does not match any variable in x, no link is added \n",
              "var1: ",var1,"\n")
    }
  }
  
  ####
  if(is.na(var2)){
    
    intercept(x) <- as.formula(paste0("~", var1))
    
  }else{
    
    if(var1 == var2){
      if(warnings){
        warning("addLink.lvm: var1 equals var2, no link is added \n",
                "var1/2: ",var1,"\n")
      }
    }
    
    
    if(var2 %in% allVars == FALSE){
      if(warnings){
        warning("addLink.lvm: var2 does not match any variable in x, no link is added \n",
                "var2: ",var2,"\n")
      }
    }
    if(var1 %in% endogenous(x) && var2 %in% endogenous(x)){
      if(missing(covariance)){
        covariance <- TRUE
      }else if(covariance == FALSE){
        if(warnings){ warning("addLink.lvm: set covariance argument to TRUE to add a covariance link to the lvm \n") }
        return(x)
      }
    }
    
    test.1 <- var1 %in% exogenous(x)
    test.2 <- var2 %in% exogenous(x)
    
    if(test.1 && test.2){
      if(warnings){
        warning("addLink.lvm: both variable are exogenous, no link is added \n",
                "var1: ",var1,"\n",
                "var2: ",var2,"\n")
      }
    }else if(test.1){
      regression(x) <- as.formula(paste(var2, var1,  sep = "~"))
    }else if(test.2){
      regression(x) <- as.formula(paste(var1, var2, sep = "~"))
    }else if(covariance){
      covariance(x) <- as.formula(paste(var1, var2, sep = "~"))  
    }else {
      if(var1 %in% endogenous(x)){
        regression(x) <- as.formula(paste(var1, var2, sep = "~"))  
      }else if(var2 %in% endogenous(x)){
        regression(x) <- as.formula(paste(var2, var1, sep = "~"))  
      }else{
        stop("unknow configuration \n")
      }
      
    }
  }
  
  return(x)
}

#' @rdname addLink
addLink.lvm.reduced <- function(x, ...){
  return(addLink.lvm(x, allVars = vars(x, lp = FALSE, xlp = TRUE) , ...))
}
# }}}


# {{{ setLink
#' @title Affect a given value to a link between two variables in a lvm
#' @name setLink
#' @description Generic interface to set a value to a link in a lvm.
#' 
#' @param x a lvm model
#' @param var1 the first variable (character) or a formula describing the link
#' @param var2 the second variable (character). Only used if var1 is a character.
#' @param value the value to which the link should be set
#' @param warnings should a warning be displayed when the link is not found in the lvm.
#' 
#' @examples
#' \dontrun{
#' setLink <- lava.penalty:::setLink
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1 ~ y2
#' 
#' m1 <- setLink(m, y3 ~ u, value = 1)
#' estimate(m1, sim(m,1e2))
#' # m1 <- setLink(m, u ~ y3, value = 1)
#' 
#' m2 <- setLink(m, y1 ~ y2, value = 0.5)
#' estimate(m2, sim(m,1e2))
#' }
`setLink` <-
  function(x,...) UseMethod("setLink")

#' @rdname setLink
setLink.lvm <- function(x, var1, var2, value, warnings = FALSE){

  res <- initVar_link(var1, var2)
  var1 <- res$var1
  var2 <- res$var2
  
  #### set the link
  if(is.na(var2)){
    intercept(x, as.formula(paste0("~",var1))) <- value
  }else if(paste(var1, var2, sep = "~") %in% coef(x)){
    regression(x, as.formula(paste(var1,var2, sep = "~"))) <- value
  }else if(paste(var1,var2, sep = ",") %in% coef(x)){
    covariance(x, as.formula(paste(var1,var2, sep = "~"))) <- value
  }else if(paste(var2,var1, sep = ",") %in% coef(x)){
    covariance(x, as.formula(paste(var1,var2, sep = "~"))) <- value
  }else{
    if(warnings){
      warning("setLink.lvm: no link was found from var1 to var2, no link is set \n",
              "var1: ",var1,"\n",
              "var2: ",var2,"\n")
    }
  }
  
  return(x)
}
# }}}


# {{{ cancel.plvm
#' @title Remove a new link between two variables in a plvm model
#' @name cancel
#' @description Generic interface to remove a link in a plvm.
#' 
#' @param x a lvm model
#' @param value the names of the links that should be removed
#' @param clean should the lvm object be simplified using the \code{clean} function
#' @param ... other argument to be passed to \code{clean}.
#' 
#' @examples 
#' 
#' #### penalized lvm ###
#' m <- lvm(Y ~ X1 + X2 + X3)
#' pm <- penalize(m)
#' cancel(pm, Y~X1+X2)
#' cancel(pm, Y~X1+X2+X3)
#'
cancel.plvm <- function(x, value, ...){

    ## normalize input
    value <- initVar_link(value, format = "txt.formula")

    ## remove penalties associated to the link
    if(any(value %in% penalty(x, type = "link"))){
        cancelPenalty(x) <- value[value %in% penalty(x, type = "link")]
    }

    ## call the method for the other classes
    x <- callS3methodParent(x, FUN = "cancel", class = "plvm", value = value, ...)
    
    return(x)
}
# }}}

