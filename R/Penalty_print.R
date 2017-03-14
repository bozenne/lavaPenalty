
# {{{ print.penaltyL12
`print.penaltyL12` <- function(x, ...){

    fctExtract <- function(x, char){
        lambda <- penalty(x, type = char)[[char]]
        if(length(lambda)>0){
            if(length(unique(lambda))==1){
                lambda <- lambda[1]
            }else{
                lambda <- "multiple values"
            }
        }else{
            lambda <- "undefined"
        }
        return(lambda)
    }

    lambda1 <- fctExtract(x, "lambda1")
    lambda2 <- fctExtract(x, "lambda2")
    lambdaG <- fctExtract(x, "lambdaG")

    ## elastic net penalty
    x.penalty <- penalty(x)

    if(x.penalty[penalty %in% c("lasso","ridge"),.N]>0){
        index.lasso <- x.penalty[penalty %in% "lasso",link]
        index.ridge <- x.penalty[penalty %in% "ridge",link]

        index.elasticNet <- intersect(index.lasso,index.ridge)
        index.lasso <- setdiff(index.lasso,index.elasticNet)
        index.ridge <- setdiff(index.ridge,index.elasticNet)
        
        if(length(index.elasticNet)>0){
            cat("Penalty: elastic net (lambda1 = ",lambda1," lambda2 = ",lambda2,") \n",
                "Links  : ", paste(index.elasticNet, collapse = " "),"\n\n",sep="")
        }
        if(length(index.lasso)>0){
            cat("Penalty: lasso (lambda1 = ",lambda1,") \n",
                "Links  : ", paste(index.lasso, collapse = " "),"\n\n",sep="")
        }
        if(length(index.ridge)>0){
            cat("Penalty: ridge (lambda2 = ",lambda2,") \n",
                "Links  : ", paste(index.ridge, collapse = " "),"\n\n",sep="")
        }                    
    }

    ## group penalty
    if(x.penalty[penalty %in% c("group lasso"),.N]>0){
        index.groupLasso <- x.penalty[penalty %in% "group lasso",link]
        group.groupLasso <- 
        all.groups <- x.penalty[penalty %in% "group lasso",unique(group)]
        
        cat("Penalty: group lasso (lambdaG = ",lambdaG,")\n",sep="")        
        sapply(all.groups, function(g){
            glink <- x.penalty[penalty=="group lasso" & group == g,link]
            cat("group ",g,": ",paste(glink, collapse = " "),"\n",sep="")
        })
        cat("\n")
    }      
}
# }}}

# {{{ print.penaltyNuclear
`print.penaltyNuclear` <- function(x, ...){
    test.nuclear <- !is.null(penalty(x, type = "link")$link)
    
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
