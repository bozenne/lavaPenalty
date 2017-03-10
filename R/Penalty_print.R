
# {{{ print.penaltyL12
`print.penaltyL12` <- function(x, ...){
    lambda1 <- penalty(x, type = "lambda1")
    if(!is.null(lambda1)){
        lambda1.mean <- mean(lambda1[lambda1>0])
    }else{
        lambda1.mean <- NA
    }

    lambda2 <- penalty(x, type = "lambda2")
    if(!is.null(lambda2)){
        lambda2.mean <- mean(lambda2[lambda2>0])
    }else{
        lambda2.mean <- NA
    }
    
    lambdaG <- penalty(x, type = "lambdaG")
    if(!is.null(lambdaG)){
        lambdaG.mean <- mean(lambdaG[lambdaG>0])
    }else{
        lambdaG.mean <- NA
    }
    ## elastic net penalty
    penalty.elasticNet <- penalty(x, type = "link", no.group = TRUE)

    if(length(penalty.elasticNet)>0){
        if(!is.na(lambda1.mean) && !is.na(lambda2.mean)){
            display.elasticNet <- paste0("elastic net (lambda1 = ",lambda1.mean," lambda2 = ",lambda2.mean,")")
        }else if(!is.na(lambda1.mean)){
            display.elasticNet <- paste0("lasso (lambda1 = ",lambda1.mean,")")
        }else if(!is.na(lambda2.mean)){
            display.elasticNet <- paste0("ridge (lambda2 = ",lambda2.mean,")")
        }else{
            display.elasticNet <- paste0("not specified")
        }        
        cat("Penalty : ", display.elasticNet,"\n",
            "Links   : ", paste(penalty.elasticNet, collapse = " "),"\n\n",sep="")
    }

    ## group penalty
    if(!is.null(penalty(x, type = "Vgroup"))){ 
        cat("Penalty : group lasso ",if(!is.na(lambdaG.mean)){paste0("(lambdaG = ",lambdaG.mean,")")},"\n")
        n.groups <- NCOL(penalty(x, type = "Vgroup"))
        sapply(1:n.groups, function(g){
            cat("group ",g,": ",paste(penalty(x, type = "link", no.elasticNet = TRUE, group = g), collapse = " "),"\n",sep="")
        })
        cat("\n")
    }      
}
# }}}

# {{{ print.penaltyNuclear
`print.penaltyNuclear` <- function(x, ...){
    test.nuclear <- !is.null(penalty(x, type = "link"))

    
    if(test.nuclear){
        allNames <- penalty(x, type = "name.reduce")
        allEndo <- penalty(x, type = "endogeneous")
        lambdaN <- penalty(x, type = "lambdaN")
        n.penaltyNuclear <- length(allNames)
        
        cat("Penalty: nuclear norm ",if(!is.null(lambdaN)){paste0("(lambdaN = ",mean(lambdaN),")")},"\n",
            "on     : ", allNames[1]," (endogenous: ",allEndo[1],")\n",sep = "")
        if(n.penaltyNuclear>1){
            sapply(1:n.penaltyNuclear, function(p){
                cat("     : ", allNames[p]," (endogeneous: ",allEndo[p],")\n",sep = "")
            })
        }
    }
    
}
# }}}

# {{{ print.plvm
#' @title Display the content of a plvm object
#
#' @param x a plvm object
#'
#' @export
`print.plvm` <- function(x, ...) {
  
    ## normal display
    out <- capture.output(lava:::print.lvm(x))
    sapply(out, function(o){cat(o,"\n")})

    ## penalty
    print(x$penalty)
    print(x$penaltyNuclear)
    
    ## export
    invisible(x)
}
# }}}

# {{{ print.plvmfit
#' @title Display the content of a plvmfit object
#
#' @param x a plvmfit object
#'
#' @export
`print.plvmfit` <- function(x,level=2,labels=FALSE,
                            coef = "penalized", lambda = "abs", only.breakpoints = TRUE, 
                            ...) {
  
  if(is.null(x$regularizationPath)){
    
    Mtempo <- CoefMat(x,labels=labels,level=level,...) 
    ncol.M <- ncol(Mtempo)

    print(Mtempo,quote=FALSE,right=TRUE)
    print(x$penalty)

    minSV <- attr(vcov(x),"minSV")
    if (!is.null(minSV) && minSV<1e-12) {
      warning("Small singular value: ", format(minSV))
    }
    pseudo <- attr(vcov(x),"pseudo")
    if (!is.null(pseudo) && pseudo) warning("Singular covariance matrix. Pseudo-inverse used.")
    
  }else if(is.null(x$regularizationPath$optimum)){
    cat("Regularization path: \n")
    print(x$regularizationPath, coef = coef, lambda = lambda, only.breakpoints = only.breakpoints)
    cat("estimated using EPSODE algorithm \n")
    
  }else{
    lava:::print.lvmfit(x)
    cat("\n Model selected using ",x$regularizationPath$optimum$criterion," \n")
    if(lambda == "abs"){
      cat("   range of lambda1.abs: ",paste(range(x$regularizationPath$performance$lambda1.abs), collapse = " "),"\n")
      cat("   best lambda1.abs    : ",x$regularizationPath$optimum$lambda1.abs,"\n")
    }else if(lambda == "nabs"){
      cat("   range of lambda1: ",paste(range(x$regularizationPath$performance$lambda1), collapse = " "),"\n")
      cat("   best lambda1    : ",x$regularizationPath$optimum$lambda1,"\n")
    }
  }
  
  invisible(x)
}

`print.regPath` <- function(x, coef = "penalized", lambda = "abs", only.breakpoints = TRUE) {
  
    printPath <- getPath(x, only.breakpoints = only.breakpoints, coefficient = coef, lambda = lambda)
    print(printPath)
    diffRow <- nrow(getPath(x)) - nrow(printPath)
    if(diffRow>0){cat("[ omitted ",diffRow," rows ] \n",sep = "")}
    
}
# }}}
