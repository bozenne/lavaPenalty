### test-estimateHook.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 10 2017 (16:42) 
## Version: 
## last-updated: feb 10 2017 (17:29) 
##           By: Brice Ozenne
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: Check the correct initialization of pLVM
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:



 m <- lvm()
 m <- regression(m,y='y1',x='x'%++%1:2)
 m <- regression(m,y='y1',x='z'%++%1:5)
 
 # simul
 set.seed(10)
 d <- sim(m,5e2)
 d <- as.data.frame(scale(d))
 
 suppressWarnings(
 start <- coef(estimate(m, data = d, control = list(iter.max = 0)))
 )
 
 # penalized lvm 
 mp <- penalize(m)
 system.time(
 lassoLVM <- estimate(mp, data = d, lambda1 = 10, control = list(trace = 2, start = start[coef(mp)], iter.max = 0), estimator = "gaussian2")
 )
 
 
 # reduced penalized lvm
 mp.red <- reduce(mp)
 system.time(
 lassoRLVM <- estimate(mp.red, data = d, lambda1 = 10, control = list(trace = 2, start = start[coef(mp.red)], iter.max = 1))
 )
 
coef(lassoRLVM)-coef(lassoLVM)[names(coef(lassoRLVM))]
lassoRLVM$opt$iterations
 



#----------------------------------------------------------------------
### test-estimateHook.R ends here
