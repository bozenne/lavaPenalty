`findNewLink` <-
  function(x,...) UseMethod("findNewLink")

`addLink` <-
  function(x,...) UseMethod("addLink")

`setLink` <-
  function(x,...) UseMethod("setLink")

`rmLink` <-
  function(x,...) UseMethod("rmLink")


#' @title Find the new possible links between variables (copied from lava::modelsearch)
#' 
#' @param x a lvm model
#' @param rm.exoexo ignore links between exogeneous variables
#' @param output return the names of the variables to link ("names") or their position ("index")
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' 
#' lava.penalty:::findNewLink(m, rm.exoexo = FALSE)
#' lava.penalty:::findNewLink(m, rm.exoexo = TRUE)
#' lava.penalty:::findNewLink(m, rm.exoexo = TRUE, output = "index")
#' 
findNewLink.lvm <- function(x, rm.exoexo, output = "names"){
  
  match.arg(output, choices = c("names","index"))
  
  AP <- with(index(x), A + t(A) + P)
  
  restricted <- c()
  for (i in seq_len(ncol(AP) - 1)){
    for (j in seq(i + 1, nrow(AP))){
      test.exo <- (rownames(AP)[i] %in% exogenous(x)) + (colnames(AP)[j] %in% exogenous(x))
      if (AP[j, i] == 0 && (rm.exoexo == FALSE || test.exo!=2)){
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

#' @title Add a new link between two variables in a lvm
#' @description Generic interface to add links to a lvm.
#' 
#' @param x a lvm model
#' @param var1 the first variable (character) or a formula describing the link to be added to the lvm
#' @param var2 the second variable (character). Only used if var1 is a character.
#' @param covariance does the link is bidirectional. Ignored if one of the variable is exogenous.
#' @param warnings should a warning be displayed when no link is added
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' 
#' lava.penalty:::addLink(m, x1 ~ y1)
#' lava.penalty:::addLink(m, y1 ~ x1)
#' coef(lava.penalty:::addLink(m, y1 ~ y2, covariance = TRUE))
#' 
#' lava.penalty:::addLink(m, "x1", "y1")
#' lava.penalty:::addLink(m, "y1", "x1")
#' coef(lava.penalty:::addLink(m, "y1", "y2", covariance = TRUE))
#' 
addLink.lvm <- function(x, var1, var2 = NA, covariance, warnings = FALSE){
 
  res <- initVar_link(var1, var2)
  var1 <- res$var1
  var2 <- res$var2
  
  if(var1 %in% vars(x) == FALSE){
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
    
    
    if(var2 %in% vars(x) == FALSE){
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
      if(var1 %in% endogenous(x) && var2 %in% latent(x)){
        regression(x) <- as.formula(paste(var1, var2, sep = "~"))  
      }else if(var2 %in% endogenous(x) && var1 %in% latent(x)){
        regression(x) <- as.formula(paste(var2, var1, sep = "~"))  
      }else{
        stop("addLink.lvm: unknow configuration \n")
      }
      
    }
  }
  
  return(x)
}

#' @title Affect a given value to a link between two variables in a lvm 
#' @description Generic interface to set a value to a link in a lvm.
#' 
#' @param x a lvm model
#' @param var1 the first variable (character) or a formula describing the link
#' @param var2 the second variable (character). Only used if var1 is a character.
#' @param value the value to which the link should be set
#' @param warnings should a warning be displayed when the link is not found in the lvm.
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1 ~ y2
#' 
#' m1 <- lava.penalty:::setLink(m, y3 ~ u, value = 1)
#' estimate(m1, sim(m,1e2))
#' # m1 <- lava.penalty:::setLink(m, u ~ y3, value = 1)
#' 
#' m2 <- lava.penalty:::setLink(m, y1 ~ y2, value = 0.5)
#' estimate(m2, sim(m,1e2))
#' 
setLink.lvm <- function(x, var1, var2 = NA, value, warnings = FALSE){
  
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

#' @title Remove a new link between two variables in a lvm model
#' @description Generic interface to remove a link in a lvm.
#' 
#' @param x a lvm model
#' @param var1 the first variable (character) or a formula describing the link
#' @param var2 the second variable (character). Only used if var1 is a character.
#' @param warnings should a warning be displayed when the link is not found in the lvm.
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1 ~ y2
#' 
#' lava.penalty:::rmLink(m, y3 ~ u)
#' 
#' coef(lava.penalty:::rmLink(m, y1 ~ y2))
rmLink.lvm <- function(x, var1, var2 = NA, warnings = FALSE){
  
  res <- initVar_link(var1, var2)
  var1 <- res$var1
  var2 <- res$var2
  
  #### remove the link
  if(is.na(var2)){
    cancel(x) <- as.formula(paste0("~",var1))
  }else if(paste(var1, var2, sep = "~") %in% coef(x)){
    cancel(x) <- as.formula(paste(var1,var2, sep = "~"))
  }else if(paste(var1,var2, sep = ",") %in% coef(x)){
    cancel(x) <-  as.formula(paste(var1,var2, sep = "~"))
  }else if(paste(var2,var1, sep = ",") %in% coef(x)){
    cancel(x) <-  as.formula(paste(var2,var1, sep = "~"))
  }else{
    if(warnings){
      warning("addLink.lvm: no link was found from var1 to var2, no link is removed \n",
              "var1: ",var1,"\n",
              "var2: ",var2,"\n")
    }
  }
  
  
  #### if unused variable remove it from the model
  if(length(grep(paste0("~",var1,"$|^",var1,"~|^",var1,"$"), x = coef(x), fixed = FALSE))==0){
      kill(x) <- as.formula(paste0(var1,"~1"))
  }
  if(!is.na(var2) && length(grep(paste0("~",var2,"$|^",var2,"~|^",var2,"$"), x = coef(x), fixed = FALSE))==0){
    kill(x) <- as.formula(paste0(var2,"~1"))
  }
  
  return(x)
}

#' @title Normalize var1 and var2
#' @description Convert var1 and var2 from formula or covariance to character
#' 
#' @examples
#' lava.penalty:::initVar_link(var1 = a~b, var2 = NA)
#' lava.penalty:::initVar_link(var1 = a ~ b, var2 = NA)
#' 
#' lava.penalty:::initVar_link(var1 = "a,b", var2 = NA)
#' lava.penalty:::initVar_link(var1 = "a", var2 = "b")
#' 
initVar_link <- function(var1, var2){
  
  if(is.na(var2) && is.character(var1)){
    if(grepl(",",var1)==TRUE){var1 <- gsub(",","~", x = var1)}
    if(grepl("~",var1)==TRUE){var1 <- as.formula(var1)}
  }
  
  if(class(var1) == "formula"){
    var2 <- all.vars(var1)[2]
    var1 <- all.vars(var1)[1]
  }
  
  return(list(var1 = var1,
              var2 = var2))
}
