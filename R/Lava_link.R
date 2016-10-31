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
#' @param rm.latent ignore links between latent variables
#' @param rm.endo ignore links between endogenous variables
#' @param output return the names of the variables to link ("names") or their position ("index")
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' lava.penalty:::findNewLink(m, rm.endo = FALSE)
#' lava.penalty:::findNewLink(m, rm.endo = TRUE)
#' 
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' 
#' lava.penalty:::findNewLink(m, rm.endo = FALSE)
#' lava.penalty:::findNewLink(m, rm.endo = TRUE)
#' lava.penalty:::findNewLink(m, rm.endo = TRUE, output = "index")
#' 
findNewLink.lvm <- function(x, rm.latent = FALSE, rm.endo = FALSE, output = "names"){
  
  match.arg(output, choices = c("names","index"))
  
  AP <- with(index(x), A + t(A) + P)
  
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

addLink.lvm.reduced <- function(x, ...){
  return(addLink.lvm(x, allVars = vars(x, lp = FALSE, xlp = TRUE) , ...))
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

#' @title Remove a new link between two variables in a lvm model
#' @description Generic interface to remove a link in a lvm.
#' 
#' @param x a lvm model
#' @param var1 the first variable (character) or a formula describing the link
#' @param var2 the second variable (character). Only used if var1 is a character.
#' @param warnings should a warning be displayed when the link is not found in the lvm.
#' @param other parameters to be passed to cancelPenalty and cancelLP
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u+x1
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1 ~ y2+y3
#' 
#' lava.penalty:::rmLink(m, y3 ~ u)
#' lava.penalty:::rmLink(m, u ~ x1)
#' lava.penalty:::rmLink(m, y1 ~ x1+u)
#' lava.penalty:::rmLink(m, u ~ x1+x2)
#' 
#' coef(lava.penalty:::rmLink(m, y1 ~ y2))
#' coef(lava.penalty:::rmLink(m, y1 ~ y2+y3))
#' 
#' # external parameter
#' parameter() <- y1 ~ y2
#' 
#' 
rmLink.lvm <- function(x, var1, var2, warnings = FALSE, ...){
  
  res <- initVar_link(var1, var2)
  var1 <- res$var1
  var2 <- res$var2
  
  #### remove the link
  if(length(var2)==0){
    cancel(x) <- as.formula(paste0("~",var1))
  }else{
    
    index.coef <- which(paste(var1, var2, sep = "~") %in% setdiff(coef(x), parameter(x)))
    index.ext <- which(paste(var1, var2, sep = "~") %in% parameter(x))
    index.cov1 <- which(paste(var1, var2, sep = ",") %in% coef(x))
    index.cov2 <- which(paste(var2, var1, sep = ",") %in% coef(x))
    
    if(length(index.coef)>0){
      f <- paste(var1,paste(var2[index.coef], collapse = " + "), sep = "~")
      cancel(x) <- as.formula(f)
    }
    if(length(index.ext)>0){
      f <- paste(var1,paste(var2[index.ext], collapse = " + "), sep = "~")
      parameter(x, remove = TRUE) <- f
    }
    if(length(index.cov1)>0){
      f <- paste(var1,paste(var2[index.cov1], collapse = " + "), sep = "~")
      cancel(x) <- as.formula(f)
    }
    if(length(index.cov2)>0){
      f <- paste(var1,paste(var2[index.cov2], collapse = " + "), sep = "~")
      cancel(x) <- as.formula(f)
    }
    
    if(length(c(index.coef,index.ext, index.cov1, index.cov2))!= length(var2) && warnings){
        warning("addLink.lvm: no link was found from var1 to var2, no link is removed \n",
                "var1: ",var1,"\n",
                "var2: ",paste(setdiff(var2, c(index.coef,index.ext, index.cov1, index.cov2)), collapse = " "),"\n")
    }
    
  }
  
  #### if unused variable remove it from the model
  if(length(var2)>0){
    length.var2 <- sapply(var2, function(var){length(grep(paste0("~",var,"$|^",var,"~|^",var,"$"), x = coef(x), fixed = FALSE))})
    kill(x) <- as.formula(paste0(paste(var2[length.var2==0],collapse = "+"),"~1"))
  }
  
  #### if penalised remove the penalty from the model
  if("plvm" %in% class(x)){
    f <- c(paste(var1,var2[index.coef], sep = "~"),
           paste(var1,var2[index.ext], sep = "~"),
           paste(var1,var2[index.cov1], sep = ","),
           paste(var1,var2[index.cov2], sep = ",")
    )
    if(any(f %in% penalty(x, type = "link"))){
      cancelPenalty(x, ...) <- f[f %in% penalty(x, type = "link")]
    }
  }
  #### if belong to a LP remove the variable from the lp
  if("lvm.reduced" %in% class(x)){
    f <- c(paste(var1,var2[index.coef], sep = "~"),
           paste(var1,var2[index.ext], sep = "~")
    )
    if(any(f %in% lp(x, type = "link"))){
      cancelLP(x, ...) <- f[f %in% lp(x, type = "link")]
    }
  }
 
  return(x)
}


# rmLink.lvm.reduced <- function(x, var1, var2, warnings = FALSE, simplify = TRUE){
#  
#   x <- rmLink.lvm(x, var1, var2, warnings = FALSE)
#   browser()
#   ## update linear predictor
#   f <- initVar_link(var1,var2, format = "txt.formula")
#   if(f %in% lp(x, type = "link")){
#    x <- cancelLP(x, link = f, simplify = simplify)  
#   }
#   
#   return(x)
# }