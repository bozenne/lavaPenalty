# NOTE: the user should not be able to remove the intercept from a model that have been reduced
`procdata` <- function(x,...) UseMethod("procdata") # to use procdata.lvm.reduced instead of procdata.lvm
deriv.lvm <- lava:::deriv.lvm # to use deriv.lvm without using :::






