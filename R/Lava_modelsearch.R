#' @export
`extendModel` <-
  function(x,...) UseMethod("extendModel")

#' @title Automatic extension of the lvm
#' 
#' @param x a lvm model
#' @param type should all links be added to the latent variable model ("all") or only thoses relevant according to a score ("modelsearch") or a LR ("modelsearchLR") test. 
#' @param covariance should covariance links be considered?
#' @param data the dataset used to perform the model search
#' @param alpha the significance threshold for retaining a new link
#' @param method the method used to adjust the p.values for multiple comparisons in the modelsearch
#' @param display.warnings should warnings be display? May occur when dealing with categorical variables or when fitting an extended model.
#' @param trace should the execution of the modelsearch be traced?
#' @param ... additional arguments to be passed to the \code{estimate} method.
#' 
#' @return a latent variable model
#' 
#' @examples 
#' mSim <- lvm()
#' regression(mSim) <- c(y1,y2,y3)~u
#' regression(mSim) <- u~x1+x2
#' categorical(mSim,labels=c("A","B","C")) <- "x2"
#' latent(mSim) <- ~u
#' covariance(mSim) <- y1~y2
#' df <- sim(mSim, 1e2)
#' 
#' m <- lvm(c(y1,y2,y3)~u)
#' latent(m) <- ~u
#' addvar(m) <- ~x1+x2 
#' 
#' coef(extendModel(m, type = "all"))
#' coef(extendModel(m, type = "all", covariance = FALSE, display.warnings = FALSE))
#' 
#' coef(extendModel(m, type = "modelsearch", data = df))
#' coef(extendModel(m, type = "modelsearchLR", data = df))
#' @export
extendModel.lvm <- function(x, type, covariance = TRUE, 
                            data, alpha = 0.05, method = "holm",
                            display.warnings = TRUE, trace = TRUE, ...){
  
  match.arg(type, choices = c("all","modelsearch", "modelsearchLR"))
  
  if(type != "all"){
   
    ## handling factor variables
    test.factor <- unlist(lapply(data, function(x){is.factor(x) + is.character(x) > 0}))
    varFactor <- names(test.factor)[test.factor]
    ls.levels <- lapply(varFactor, function(x){levels(data[[x]])})
    names(ls.levels) <- varFactor
    
    lvmfit <- estimate(x, data = data, ...)
    cv <- FALSE
    
    while(cv == FALSE){ 
      if(trace){cat("*")}
      
      if(display.warnings){
        resSearch <- do.call(type, args = list(lvmfit, na.omit = TRUE, silent = TRUE)) 
      }else{
        suppressWarnings(
          resSearch <- do.call(type, args = list(lvmfit, na.omit = TRUE, silent = TRUE))
        )
      }
      
      if(covariance){
        index <- 1:length(resSearch$var)
      }else{
        index <- which(unlist(lapply(resSearch$var, function(var){sum(var %in% endogenous(x))}))<2)
      }
      
      if(tail(p.adjust(resSearch$test[index,"P-value"], method = method), 1) < alpha){
        var1 <- tail(resSearch$var[index],1)[[1]][,1]
        var2 <- tail(resSearch$var[index],1)[[1]][,2]
        
        if(var1 %in% vars(x) == FALSE){
          var1 <- renameFactor(var1, ls.levels = ls.levels)
          if(display.warnings && length(ls.levels[[var1]])>2){
            warning("extendModel.lvm: one of the levels of a factor variable reach the significance level \n",
                    "a link with the whole variable is added in the model \n")
          }
        }
        if(var2 %in% vars(x) == FALSE){
          var2 <- renameFactor(var2, ls.levels = ls.levels)
          if(display.warnings && length(ls.levels[[var2]])>2){
            warning("extendModel.lvm: one of the levels of a factor variable reach the significance level \n",
                    "a link with the whole variable is added in the model \n")
          }
        }
        x <- addLink(x, var1, var2,
                     covariance = covariance, warnings = display.warnings)
        
        lvmfit <- estimate(x, data = data, ...)
      }else{
        cv <- TRUE
      }
      
      if(lvmfit$opt$iterations == lvmfit$control$iter.max){
        return(NA)
      }
    }
    if(trace){cat("\n")}
    
    return(lvmfit)
    
  }else{ # all
    
    newlinks <- findNewLink(x, rm.exoexo = TRUE)
    for(iterLink in 1:nrow(newlinks)){
     x <- addLink(x, newlinks[iterLink,1], newlinks[iterLink,2], covariance = covariance,
                  warnings = FALSE)
    }
    
    return(x)
  }
  

}

#' @title Model searching using a likelihood ratio test
#' 
#' @param x a lvm model
#' @param na.omit do not export the results for links where the extended model has not converged
#' @param display.warnings should the warnings encountered after the fit of an extended model be displayed?
#' @param silent should the execution of the modelsearch be traced?
#' @param ... additional arguments to be passed to the \code{estimate} method.
#' 
#' @return an object of class lvmfit
#' 
#' @examples 
#' mSim <- lvm()
#' regression(mSim) <- c(y1,y2,y3)~u
#' regression(mSim) <- u~x1+x2
#' categorical(mSim,labels=c("A","B","C")) <- "x2"
#' latent(mSim) <- ~u
#' covariance(mSim) <- y1~y2
#' df <- sim(mSim, 1e2)
#' 
#' m <- lvm(c(y1,y2,y3)~u)
#' latent(m) <- ~u
#' addvar(m) <- ~x1+x2 
#' 
#' fit <- estimate(m, data = df)
#' modelsearchLR(fit)
#' modelsearchLR(fit, na.omit = TRUE)
#' 
#' @export

`modelsearchLR` <- function(object, ...) UseMethod("modelsearchLR")

modelsearchLR.lvmfit <- function (object, na.omit = FALSE, display.warnings = FALSE, silent = FALSE, ...){

  #### newlinks 
  restricted <- findNewLink(object$model, rm.exoexo = FALSE, output = "names")
  seq_i <- seq_len(NROW(restricted))
  
  #### initialisation
  M.test <- cbind("Test Statistic" = rep(NA,length(seq_i)),
                  "P-value" = rep(NA,length(seq_i))
  )
  ls.var <- list()
  
  if(silent == FALSE){pb <- utils::txtProgressBar(max = tail(seq_i,1), style = 3) }
  
  for (iterI in seq_i) {
    
    newmodel <- addLink(object$model, var1 = restricted[iterI,1], var2 = restricted[iterI,2],
                        covariance = all(restricted[iterI,1:2] %in% endogenous(object$model)))
    newcontrol <- object$control
    newcontrol$start <- coef(object)
    newcontrol$trace <- FALSE
    
    if(display.warnings){
      newfit <- tryCatch(estimate(newmodel, data = object$data$model.frame, control = newcontrol, 
                                  missing = "lvm.missing" %in% class(object), ...),
                         error = function(x){NA},
                         finally = function(x){x})
    }else{
      newfit <- suppressWarnings(tryCatch(estimate(newmodel, data = object$data$model.frame, control = newcontrol, 
                                                   missing = "lvm.missing" %in% class(object), ...),
                                          error = function(x){NA},
                                          finally = function(x){x}))
    }
    browser()
    
    if("lvmfit" %in% class(newfit)){ # test lvmfit is not an error
      if(newfit$opt$convergence == 0){ # test whether lvmfit has correctly converged
        compareT <- compare(object,newfit)
        M.test[iterI,] <- c(compareT$statistic[[1]], compareT$p.value[[1]])
      }
    }
    
    ls.var[[iterI]] <- matrix(c(restricted[iterI,1], restricted[iterI,2]), nrow = 1)
    if(silent == FALSE){ utils::setTxtProgressBar(pb, value = iterI) }
    
  }
  if(silent == FALSE){  close(pb) }
  
  if(na.omit){
    index.NNA <- which(rowSums(!is.na(M.test))>0)
    M.test <- M.test[index.NNA,,drop = FALSE]
    ls.var <- ls.var[index.NNA]
  }
  
  #### reorder
  index.order <- order(M.test[,"P-value"], decreasing = TRUE)
  ls.var <- ls.var[index.order]
  M.test <- M.test[index.order,]
  
  #### export
  M.res <- cbind(apply(M.test, 2, signif, digit = 4), unlist(lapply(ls.var, paste, collapse = ",")))
  colnames(M.res) <- c("Score:", "S P(S>s)", "Index")
  rownames(M.res) <- rep("",nrow(M.res))
  
  output <- list(res = M.res,
                 test = M.test,
                 var = ls.var)
  class(output) <- "modelsearch"
  return(output)
}
