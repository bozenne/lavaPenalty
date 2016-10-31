#' @title Display the content of a plvm object
#
#' @param x a plvm object
#'
#' @export
`print.plvm` <- function(x, ...) {
  
  ## normal display
  out <- capture.output(lava:::print.lvm(x))
  # if(!is.null(x$penaltyNuclear$name.Y)){
  #   browser()
  #   charY <- paste0(x$penaltyNuclear$name.Y," ~ ")
  #   indexEq <- grep(charY,out)
  #   out[indexEq] <- gsub(charY, replacement = paste0(charY,LCSseq(x$penaltyNuclear$name.X),"(image)+"),x = out[indexEq])
  # }
  sapply(out, function(o){cat(o,"\n")})
  
  ## additional display - lasso
  if(!is.null(penalty(x, type = "link"))){
    if(penalty(x, type = "lambda1")>0 && penalty(x, type = "lambda2")>0){
      penaltyType <- "Elastic net"
    }else if(penalty(x, type = "lambda1")>0){
      penaltyType <- "Lasso"
    }else if(penalty(x, type = "lambda1")>0){
      penaltyType <- "Ridge"
    }else{
      penaltyType <- "None"
    }
    
    if(all(penalty(x, type = "group")<1)){
      cat("Penalty: ", penaltyType,"\n",
          "On     : ", paste(penalty(x, type = "link"), collapse = " "),"\n")
    }else{
      
      test.lasso <- (penalty(x, type = "group")<1)
      if(any(test.lasso==1)){
        cat("Penalty: ", penaltyType,"\n",
            "On     : ", paste(penalty(x, type = "link", group = 0), collapse = " "),"\n")
      }
      
      ls.penalty <- tapply(penalty(x, type = "link")[test.lasso!=1], penalty(x, type = "group")[test.lasso!=1],list)
      cat("Penalty: Grouped lasso \n")
      lapply(ls.penalty, function(x){cat("On     :",paste(x, collapse = " "),"\n")})
    }
    cat("\n")
  }
  if(!is.null(x$penaltyNuclear$name.Y)){
    penaltyType <- "Nuclear norm"
    cat("Penalty: ", penaltyType,"\n",
        "on     : ", LCSseq(x$penaltyNuclear$name.X)," (outcome: ",x$penaltyNuclear$name.Y,")\n",sep = "")
  }
  ## export
  invisible(x)
}

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
  
    if(penalty(x, type = "lambda1")>0 || penalty(x, type = "lambda2")>0){
      Mtempo <- rbind(Mtempo, "Penalization:" = rep("", ncol.M))
    }
    if(penalty(x, type = "lambda1")>0){
      Mtempo <- rbind(Mtempo, "   L1 lambda (abs)" = c(penalty(x, type = "lambda1.abs"), rep("",ncol.M-1)))
      Mtempo <- rbind(Mtempo, "   L1 lambda" = c(penalty(x, type = "lambda1"), rep("",ncol.M-1)))
    }
    if(penalty(x, type = "lambda2")>0){
      Mtempo <- rbind(Mtempo, "   L2 lambda (abs)" = c(penalty(x, type = "lambda2.abs"), rep("",ncol.M-1)))
      Mtempo <- rbind(Mtempo, "   L2 lambda" = c(penalty(x, type = "lambda2"), rep("",ncol.M-1)))
    }
    
    print(Mtempo,quote=FALSE,right=TRUE)
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

  
