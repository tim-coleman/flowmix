### NB: set the path to wherever /flowmix_CVXR/fused_lasso_C/ is
### e.g., setwd("C:/Users/drain/Box Sync/USC_Stuff/Featureless_Flowmix/fused_lasso_C")

setwd("C:/Users/drain/Box Sync/USC_Stuff/Featureless_Flowmix/fused_lasso_C")
dyn.load("tf_dp_ryan_orig.dll")

fused_lasso = function(z, lam) {
  if (!is.loaded("prox_dp_R")) {
    dyn.load("prox_R.so")
  }
  
  o = .C("prox_dp_R",
         n=as.integer(length(z)),
         z=as.double(z),
         lam=as.double(lam),
         beta=as.double(numeric(length(z))),
         dup=FALSE)
  
  return(o$beta)
}
