### Penalty_logLik.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: mar  7 2017 (14:47) 
## Version: 
## last-updated: mar 15 2017 (18:26) 
##           By: Brice Ozenne
##     Update #: 12
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

`logLik.plvmfit` <- function(object, addPenalty = TRUE, ...){
    resLogLik <- lava.reduce::callS3methodParent(object,"logLik")
    if(addPenalty){
        resLogLik[1] <-  resLogLik + object$opt$objectivePenalty(coef(object))
    }
    attr(resLogLik, "df") <- NA    
    # [TODO] df can be easily computed for the lasso but otherwise ...
    
    return(resLogLik)
}


#----------------------------------------------------------------------
### Penalty_logLik.R ends here
